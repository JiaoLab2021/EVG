#include "../include/vcf_convert.hpp"

#include <exception>
#include <stdexcept>

#include <htslib/kseq.h>

#include "../include/hts_file.hpp"
#include "../include/logger.hpp"

using namespace std;

KSEQ_INIT(InputFile*, read_input)

// Store the conversion inputs, thresholds and thread count.
Convert::Convert(
    string refFileName,
    string vcfFileName,
    int readLen,
    double MAF,
    double MISSRATE,
    string outFileName,
    uint32_t threads
) : refFileName_(refFileName), vcfFileName_(vcfFileName), readLen_(readLen), MAF_(MAF), MISSRATE_(MISSRATE), outFileName_(outFileName), threads_(threads) {}

// Load the reference into memory, build the ##contig header lines and write CHROMOSOME.NAME (used by GraphTyper2).
void Convert::build_reference_index() {
    ofstream outFile;
    outFile.open("CHROMOSOME.NAME", ios::out);
    if(!outFile.is_open()) {
        throw runtime_error("'CHROMOSOME.NAME': No such file or directory.");
    }

    InputFile reference(refFileName_, threads_);
    kseq_t *ks = kseq_init(&reference);

    int status;
    while((status = kseq_read(ks)) >= 0) {
        string chromosome = ks->name.s;
        uint32_t chrLen = ks->seq.l;
        string sequence = ks->seq.s;

        refIndexS_.chrLenTxt += "##contig=<ID=" + chromosome + ",length=" + to_string(chrLen) + ">\n";
        refIndexS_.chrSeqMap[chromosome] = sequence;
        refIndexS_.chrLenMap[chromosome] = chrLen;

        outFile << chromosome + "\n";
    }
    if (status < -1) {
        throw runtime_error("'" + refFileName_ + "' -> truncated or damaged FASTA/FASTQ file.");
    }
    kseq_destroy(ks);

    outFile.close();
}

// Convert the VCF: workers run convert_record; this loop keeps the ordered steps (header, column count, merging, two pending records) and writes.
void Convert::convert_vcf() {
    const map<string, string>& chrSeqMap = refIndexS_.chrSeqMap;

    SAVE SAVECLASS(outFileName_, threads_);

    stringstream outStream;
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10;
    outStream.str().reserve(CACHE_SIZE);

    VCFReader<ConvertResult> VCFOPENCLASS(vcfFileName_, threads_, [this](VCFINFOSTRUCT& INFOSTRUCTTMP, ConvertResult& result) {
        convert_record(INFOSTRUCTTMP, result);
    });

    uint32_t lineNum = 0;

    ProcessedRecord<ConvertResult> record;
    ProcessedRecord<ConvertResult> prePre;
    ProcessedRecord<ConvertResult> pre;

    // Write a pending record: the line built on the worker, or the line built now if it is not ready.
    auto write_pending = [&](ProcessedRecord<ConvertResult>& pending) {
        if (pending.value.lineReady) {
            outStream << pending.value.outLine;
        } else {
            outStream << format_record(pending.info);
        }
    };

    while (VCFOPENCLASS.read(record)) {
        VCFINFOSTRUCT& INFOSTRUCTTMP = record.info;

        uint32_t thisLineNum = INFOSTRUCTTMP.lineVec.size();

        // ------ Header: ##contig lines before #CHROM, old length lines dropped (BayesTyper) ------
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            if (INFOSTRUCTTMP.line.find("#CHROM") != string::npos) {
                outStream << refIndexS_.chrLenTxt;

                outStream << INFOSTRUCTTMP.line + "\n";

                lineNum = thisLineNum;
            }
            else if (INFOSTRUCTTMP.line.find(",length") == string::npos) {
                outStream << INFOSTRUCTTMP.line + "\n";
            }

            continue;
        }

        // ------ Column count; the worker's messages and error count only for records that pass it ------
        if (lineNum > 0 && lineNum != thisLineNum) {
            warning_stream() << "The number of columns in the VCF file is inconsistent. Skipping this site -> "
                 << INFOSTRUCTTMP.CHROM << " "
                 << record.value.pos << endl;
            continue;
        }

        // ------ The worker's messages, its error, then its verdict ------
        if (!record.get().keep) {
            continue;
        }

        // ------ Merge records at one position; two pending records restore the order after shifting ------
        if (INFOSTRUCTTMP.CHROM == pre.info.CHROM && INFOSTRUCTTMP.POS == pre.info.POS) {
            merge_loci(pre.info, INFOSTRUCTTMP, chrSeqMap);
            pre.value.lineReady = false;
            continue;
        } else if (INFOSTRUCTTMP.CHROM == prePre.info.CHROM && INFOSTRUCTTMP.POS == prePre.info.POS) {
            merge_loci(prePre.info, INFOSTRUCTTMP, chrSeqMap);
            prePre.value.lineReady = false;
            continue;
        } else if (INFOSTRUCTTMP.CHROM == prePre.info.CHROM && INFOSTRUCTTMP.POS > prePre.info.POS && INFOSTRUCTTMP.POS < pre.info.POS) {
            record.swap(pre);
        } else if (INFOSTRUCTTMP.CHROM == prePre.info.CHROM && INFOSTRUCTTMP.POS < prePre.info.POS) {
            warning_stream() << "The variant's start position is less than the previous one. Skipping this site -> "
                << INFOSTRUCTTMP.CHROM << " "
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // ------ Write the oldest pending record ------
        if (INFOSTRUCTTMP.CHROM != prePre.info.CHROM || prePre.info.POS != INFOSTRUCTTMP.POS) {
            if (!prePre.info.CHROM.empty()) {
                write_pending(prePre);
            }
            prePre.swap(pre);
            pre.swap(record);
        }

        if (outStream.tellp() >= CACHE_SIZE) {
            string outTxt = outStream.str();
            SAVECLASS.save(outTxt);
            outStream.str(string());
            outStream.clear();
        }
    }

    // ------ Write the last two pending records ------
    write_pending(prePre);
    write_pending(pre);

    if (outStream.tellp() > 0) {
        string outTxt = outStream.str();
        SAVECLASS.save(outTxt);
        outStream.str(string());
        outStream.clear();
    }
}

// Run the per-record steps (GT, SNP filter, reference checks, INFO, alleles, genotypes) on a worker; set keep and build the output line early.
void Convert::convert_record(VCFINFOSTRUCT& INFOSTRUCTTMP, ConvertResult& result) const {
    result.pos = INFOSTRUCTTMP.POS;

    static const std::regex endReg(R"(END=\d+)");
    static const std::regex svlenReg1(R"(SVLEN=[-\d+,]*;)");
    static const std::regex svlenReg2(R"(;SVLEN=[-\d+,]*)");
    static const std::regex svlenReg3(R"(SVLEN=[-\d+,]*)");

    const map<string, string>& chrSeqMap = refIndexS_.chrSeqMap;
    const map<string, uint32_t>& chrLenMap = refIndexS_.chrLenMap;

    uint32_t thisLineNum = INFOSTRUCTTMP.lineVec.size();

    // ------ GT field ------
    vector<string> formatVec = split(INFOSTRUCTTMP.FORMAT, ":");
    vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");
    uint32_t gtIndex = distance(formatVec.begin(), gtItera);
    if (gtIndex == formatVec.size()) {
        warning_stream() << "The FORMAT field does not contain genotype (GT) information at the following location: "
             << INFOSTRUCTTMP.CHROM << ":" << INFOSTRUCTTMP.POS << endl;
        return;
    }

    // ------ SNP filter by MAF and missing rate, only when both thresholds are set ------
    if (MAF_ > 0 && MISSRATE_ < 1) {
        INFOSTRUCTTMP.ID = VCFOPEN::get_TYPE(
            INFOSTRUCTTMP.LEN,
            INFOSTRUCTTMP.ALTVec
        );

        if (INFOSTRUCTTMP.ID == "SNP") {
            double MAFTMP;
            double MISSRATETMP;

            unordered_map<int, vector<string> > GTVecMapTmp;
            uint32_t misSampleNum;
            tie(GTVecMapTmp, misSampleNum) = VCFOPEN::get_gt(INFOSTRUCTTMP.lineVec);

            if (GTVecMapTmp.size() - misSampleNum <= 1) return;

            tie(MAFTMP, MISSRATETMP) = VCFOPEN::calculate(
                INFOSTRUCTTMP.lineVec,
                misSampleNum
            );

            if (MAFTMP < MAF_ || MISSRATETMP > MISSRATE_) return;
        }
    }

    // ------ Chromosome, symbolic ALT, read length (Paragraph) and REF checks ------
    if (chrSeqMap.find(INFOSTRUCTTMP.CHROM) == chrSeqMap.end() || chrLenMap.find(INFOSTRUCTTMP.CHROM) == chrLenMap.end()) {
        warning_stream() << "The chromosome '"
            << INFOSTRUCTTMP.CHROM
            << "' is not found in the reference file '" << refFileName_ << "'"
            << endl;
        return;
    }

    if (INFOSTRUCTTMP.lineVec[4].find("<") != string::npos) {
        warning_stream() << "The query sequence contains the symbol '<'. Skipping this site -> "
            << INFOSTRUCTTMP.CHROM << " "
            << INFOSTRUCTTMP.POS << " "
            << INFOSTRUCTTMP.lineVec[4] << endl;
        return;
    }

    if (static_cast<uint32_t>(readLen_) > INFOSTRUCTTMP.POS) {
        warning_stream() << "The variant's start position is less than the read length. Skipping this site -> "
            << INFOSTRUCTTMP.CHROM << " "
            << INFOSTRUCTTMP.POS << endl;
        return;
    }

    if (chrLenMap.at(INFOSTRUCTTMP.CHROM) < (INFOSTRUCTTMP.END)) {
        warning_stream() << "The variant's end position exceeds the length of the chromosome -> "
            << INFOSTRUCTTMP.CHROM << " "
            << INFOSTRUCTTMP.POS << endl;
        return;
    }

    string trueRefSeq = chrSeqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
    if (trueRefSeq != INFOSTRUCTTMP.REF) {
        warning_stream() << "Sequence discrepancy between the reference genome and the VCF, replacing with the sequence from the reference genome -> "
            << INFOSTRUCTTMP.CHROM << " "
            << INFOSTRUCTTMP.POS << endl;

        INFOSTRUCTTMP.lineVec[3] = trueRefSeq;
        INFOSTRUCTTMP.REF = trueRefSeq;
    }

    // ------ INFO: END from the REF length (Paragraph), SVLEN removed (GraphTyper2) ------
    INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], endReg, "END=" + to_string(INFOSTRUCTTMP.END));

    INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], svlenReg1, "");
    INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], svlenReg2, "");
    INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], svlenReg3, ".");

    // ------ Alleles: shared first base, no REF copies, ACGTN only, no repeats (BayesTyper) ------
    process_sequences(INFOSTRUCTTMP, chrSeqMap);

    check_repeated_alleles(INFOSTRUCTTMP);

    // ------ Genotypes: GT only, phased and diploid (PanGenie) ------
    uint32_t maxGT = 0;
    INFOSTRUCTTMP.lineVec[8] = "GT";
    if (gtItera != formatVec.end()) {
        for (size_t i = 9; i < thisLineNum; i++) {
            string gt;

            string genotypeInfo = INFOSTRUCTTMP.lineVec[i];
            if (genotypeInfo.empty()) {
                warning_stream() << "Genotype (GT) is empty ->"
                    << INFOSTRUCTTMP.CHROM << " "
                    << INFOSTRUCTTMP.POS << endl;
                gt = "0|0";
            } else {
                gt = VCFOPEN::sample_gt(INFOSTRUCTTMP.lineVec, i, gtIndex);
            }

            if (gt == ".") {
                gt = ".|.";
            } else if (gt == "0") {
                gt = "0|0";
            }

            if (gt.find("/") != string::npos) {
                static const std::regex reg("/");
                gt = regex_replace(string(gt), reg, string("|"));
            }

            vector<string> gtVec = VCFOPEN::gt_split(gt);
            for (auto it1 : gtVec) {
                try {
                    if (stoul(it1) > maxGT) {
                        maxGT = stoi(it1);
                    }
                } catch(const std::invalid_argument& e) {
                    continue;
                }
            }

            if (gtVec.size() == 1) {
                gt = gtVec[0] + "|0";
            } else if (gtVec.size() > 2) {
                gt = gtVec[0] + "|" + gtVec[1];
            }

            INFOSTRUCTTMP.lineVec[i] = gt;
        }
    } else {
        warning_stream() << "Genotype (GT) not found in FORMAT column. Skipping this site -> "
             << INFOSTRUCTTMP.CHROM << " "
             << INFOSTRUCTTMP.POS << endl;
        return;
    }

    // ------ Genotype allele beyond the ALT count (PanGenie), variant ending at the chromosome end (Paragraph) ------
    if (INFOSTRUCTTMP.ALTVec.size() < maxGT) {
        warning_stream() << "The number of genotypes (GTs) exceeds the number of query sequences at: "
             << INFOSTRUCTTMP.CHROM << " "
             << INFOSTRUCTTMP.POS << endl;
        return;
    }

    if (chrLenMap.at(INFOSTRUCTTMP.CHROM) == (INFOSTRUCTTMP.END)) {
        warning_stream() << "The variant's end position equals the chromosome length. Skipping this site -> "
             << INFOSTRUCTTMP.CHROM << " "
             << INFOSTRUCTTMP.POS << endl;
        return;
    }

    if (INFOSTRUCTTMP.line.empty()) return;

    result.keep = true;

    // ------ Output line if no later record merges into this one; built at write time instead if building it prints or throws ------
    thread_local VCFINFOSTRUCT trimmed;
    trimmed = INFOSTRUCTTMP;
    try {
        DiagnosticCapture quiet;
        result.outLine = format_record(trimmed);
        result.lineReady = quiet.size() == 0;
    } catch (...) {
        result.lineReady = false;
    }
    if (!result.lineReady) result.outLine.clear();
}

// Remove ALT alleles absent from the genotypes and return the output line ("" if the ALT alleles repeat or none is left).
string Convert::format_record(VCFINFOSTRUCT& INFOSTRUCTTMP) const {
    trim_update_alt_alleles(INFOSTRUCTTMP);
    check_repeated_alleles(INFOSTRUCTTMP);
    if (INFOSTRUCTTMP.line.empty() || INFOSTRUCTTMP.ALTVec.empty()) {
        return string();
    }
    return join(INFOSTRUCTTMP.lineVec, "\t") + "\n";
}

// Merge a record into the pending record at its position: extend REF, append ALTs, renumber genotypes (GT only: field 0 of each sample).
void Convert::merge_loci(
    VCFINFOSTRUCT& preINFOSTRUCTTMP,
    VCFINFOSTRUCT& INFOSTRUCTTMP,
    const std::map<std::string, std::string>& chrSeqMap
) const {
    int32_t preAltNum = preINFOSTRUCTTMP.ALTVec.size();

    bool isPreLonger = preINFOSTRUCTTMP.REF.size() > INFOSTRUCTTMP.REF.size();
    string REF = isPreLonger ? preINFOSTRUCTTMP.REF : INFOSTRUCTTMP.REF;
    string diffREF = REF.substr(std::min(preINFOSTRUCTTMP.REF.size(), INFOSTRUCTTMP.REF.size()));

    unordered_map<int, vector<string> > preGTVecMapTmp, GTVecMapTmp;
    uint32_t misSampleNum;
    tie(preGTVecMapTmp, misSampleNum) = VCFOPEN::get_gt(preINFOSTRUCTTMP.lineVec);
    tie(GTVecMapTmp, misSampleNum) = VCFOPEN::get_gt(INFOSTRUCTTMP.lineVec);

    for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
        vector<string> preGtVec = preGTVecMapTmp.at(i - 9);
        vector<string> gtVec = GTVecMapTmp.at(i - 9);

        vector<string> preGtVecTmp = sort_unique(preGtVec);
        vector<string> gtVecTmp = sort_unique(gtVec);

        if ((preGtVecTmp.back() != "0" && preGtVecTmp.back() != ".") || (gtVecTmp.back() == "0" || gtVecTmp.back() == ".")) {
            continue;
        }

        for (size_t j = 0; j < preGtVec.size(); j++) {
            if (gtVec[j] == "0" || gtVec[j] == "."){
                continue;
            }
            preGtVec[j] = to_string(stoul(gtVec[j]) + preAltNum);
        }
        vector<string> sampleVec = split(preINFOSTRUCTTMP.lineVec[i], ":");
        sampleVec[0] = join(preGtVec, "|");
        preINFOSTRUCTTMP.lineVec[i] = join(sampleVec, ":");
    }

    preINFOSTRUCTTMP.REF = REF;
    preINFOSTRUCTTMP.lineVec[3] = REF;

    preINFOSTRUCTTMP.LEN = preINFOSTRUCTTMP.REF.size();

    if (isPreLonger) {
        for (auto & altSeq : INFOSTRUCTTMP.ALTVec) {
            altSeq += diffREF;
        }
    } else {
        for (auto & altSeq : preINFOSTRUCTTMP.ALTVec) {
            altSeq += diffREF;
        }
    }
    preINFOSTRUCTTMP.ALTVec.insert(preINFOSTRUCTTMP.ALTVec.end(), INFOSTRUCTTMP.ALTVec.begin(), INFOSTRUCTTMP.ALTVec.end());
    preINFOSTRUCTTMP.ALT = join(preINFOSTRUCTTMP.ALTVec, ",");
    preINFOSTRUCTTMP.lineVec[4] = preINFOSTRUCTTMP.ALT;

    preINFOSTRUCTTMP.END = preINFOSTRUCTTMP.POS + preINFOSTRUCTTMP.REF.size() - 1;

    static const std::regex endReg(R"(END=\d+)");
    preINFOSTRUCTTMP.lineVec[7] = regex_replace(preINFOSTRUCTTMP.lineVec[7], endReg, "END=" + to_string(preINFOSTRUCTTMP.END));

    process_sequences(preINFOSTRUCTTMP, chrSeqMap);
}

// Give REF and ALT the same first base, and clear the record if an ALT equals REF or has non-ACGTN bases.
void Convert::process_sequences(VCFINFOSTRUCT& INFOSTRUCTTMP, const std::map<std::string, std::string>& chrSeqMap) const {
    for (size_t i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++) {
        string qrySeq = INFOSTRUCTTMP.ALTVec[i];

        if (qrySeq[0] != INFOSTRUCTTMP.REF[0] && (qrySeq.length() > 1 || INFOSTRUCTTMP.REF.length() > 1)) {
            INFOSTRUCTTMP.POS = INFOSTRUCTTMP.POS - 1;
            INFOSTRUCTTMP.lineVec[1] = to_string(INFOSTRUCTTMP.POS);
            INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.LEN + 1;
            INFOSTRUCTTMP.REF = chrSeqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
            INFOSTRUCTTMP.lineVec[3] = INFOSTRUCTTMP.REF;

            for (size_t j = 0; j < INFOSTRUCTTMP.ALTVec.size(); j++) {
                INFOSTRUCTTMP.ALTVec[j] = INFOSTRUCTTMP.REF[0] + INFOSTRUCTTMP.ALTVec[j];
            }

            qrySeq = INFOSTRUCTTMP.ALTVec[i];
        }

        if (qrySeq == INFOSTRUCTTMP.REF) {
            warning_stream() << "The REF and ALT sequences are the same. Skipping this site -> "
                 << INFOSTRUCTTMP.CHROM << " "
                 << INFOSTRUCTTMP.POS << endl;
            INFOSTRUCTTMP.line.clear();
            break;
        }

        smatch results;
        static const std::regex atgcReg("[^ATGCNatgcn]");
        if (regex_search(qrySeq, results, atgcReg)) {
            warning_stream() << "Sequence contains non-ATGCNatgcn characters, skip this site -> "
                 << INFOSTRUCTTMP.CHROM << " "
                 << INFOSTRUCTTMP.POS << endl;
            INFOSTRUCTTMP.line.clear();
            break;
        }
    }

    INFOSTRUCTTMP.lineVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");
}

// Remove ALT alleles absent from every genotype and renumber the genotypes (samples hold GT only: field 0).
void Convert::trim_update_alt_alleles(VCFINFOSTRUCT& INFOSTRUCTTMP) const {
    const uint32_t gtIndex = 0;

    set<uint32_t> allelesInGenotypes;
    allelesInGenotypes.insert(0);

    for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
        string gt = split(INFOSTRUCTTMP.lineVec[i], ":")[gtIndex];

        vector<string> alleles = VCFOPEN::gt_split(gt);

        for (const string& allele : alleles) {
            try {
                allelesInGenotypes.insert(stoul(allele));
            } catch (const std::invalid_argument& e) {
            }
        }
    }

    if (allelesInGenotypes.size() == INFOSTRUCTTMP.ALTVec.size() + 1) {
        return;
    }

    vector<string> newALTVec;
    set<uint32_t> allelesToRemove;
    for (size_t i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++) {
        uint32_t alleleIndex = i + 1;
        if (allelesInGenotypes.count(alleleIndex) == 0) {
            allelesToRemove.insert(alleleIndex);
        } else {
            newALTVec.push_back(INFOSTRUCTTMP.ALTVec[i]);
        }
    }

    INFOSTRUCTTMP.ALTVec = newALTVec;
    INFOSTRUCTTMP.lineVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");

    for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
        vector<string> gtVec = split(INFOSTRUCTTMP.lineVec[i], ":");
        string gt = gtVec[gtIndex];

        vector<string> alleles = VCFOPEN::gt_split(gt);

        for (string& allele : alleles) {
            try {
                uint32_t oldAlleleNumber = stoul(allele);
                uint32_t newAlleleNumber = oldAlleleNumber;

                for (const auto& alleleToRemove : allelesToRemove) {
                    if (alleleToRemove < oldAlleleNumber) {
                        --newAlleleNumber;
                    }
                }

                allele = to_string(newAlleleNumber);
            } catch (const std::invalid_argument& e) {
            }
        }

        gtVec[gtIndex] = join(alleles, "|");
        INFOSTRUCTTMP.lineVec[i] = join(gtVec, ":");
    }
}

// Clear the record if an ALT allele is repeated (BayesTyper rejects repeats).
void Convert::check_repeated_alleles(VCFINFOSTRUCT& infoStruct) const {
    for (const auto& it : infoStruct.ALTVec) {
        if (count(infoStruct.ALTVec.begin(), infoStruct.ALTVec.end(), it) > 1) {
            warning_stream() << "Repeated alleles found. Skipping this site -> "
                << infoStruct.CHROM << " "
                << infoStruct.POS << endl;
            infoStruct.line.clear();
            break;
        }
    }
}

// Return the alleles sorted numerically without duplicates.
vector<string> Convert::sort_unique(vector<string> vec) const {
    sort(vec.begin(), vec.end(), [](const string& a, const string& b) {
        return stoul(a) < stoul(b);
    });
    auto last = unique(vec.begin(), vec.end());
    vec.erase(last, vec.end());
    return vec;
}
