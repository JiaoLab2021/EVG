#include "../include/vcf_split.hpp"

#include <regex>
#include <stdexcept>

#include "../include/logger.hpp"

using namespace std;

// Store the inputs, the search distance, the FILTER switch, the output prefix and the thread count.
VCFSplit::VCFSplit(
    string vcfFileName,
    string baseVcfFileName,
    int length,
    bool filterBool,
    string prefix,
    uint32_t threads
) : vcfFileName_(vcfFileName), baseVcfFileName_(baseVcfFileName), length_(length), filterBool_(filterBool), prefix_(prefix), threads_(threads) {}

// Find the REF and ALT lengths of one data record on a worker.
void VCFSplit::measure_record(VCFINFOSTRUCT& INFOSTRUCTTMP, SplitResult& result) const
{
    result.keep = get_len(INFOSTRUCTTMP, result.refLen, result.qryLen);
}

// Get the REF and ALT lengths of a called record; return false for reference or missing genotypes, non-PASS records with -f, or unknown lengths.
bool VCFSplit::get_len(
    const VCFINFOSTRUCT& INFOSTRUCTTMP,
    uint32_t& refLen,
    uint32_t& qryLen
) const
{
    vector<string> formatVec = split(strip(INFOSTRUCTTMP.lineVec[8], '\n'), ":");
    vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");

    string gt;
    int gtIndex = 0;
    if (gtItera != formatVec.end()) {
        gtIndex = distance(formatVec.begin(), gtItera);
        gt = strip(VCFOPEN::sample_gt(INFOSTRUCTTMP.lineVec, INFOSTRUCTTMP.lineVec.size() - 1, gtIndex), '\n');
    } else {
        throw runtime_error("GT not in FORMAT column -> " + INFOSTRUCTTMP.line);
    }

    if (gt == "0/0" || gt == "." || gt == "./." || gt == ".|." || gt == "0|0") {
        return false;
    }

    if (INFOSTRUCTTMP.FILTER != "PASS" && filterBool_) {
        return false;
    }

    vector<string> evaluate_gt;
    if (filterBool_) {
        evaluate_gt = VCFOPEN::gt_split(gt);
    }
    else {
        evaluate_gt = {"1", "1"};
    }

    if (INFOSTRUCTTMP.ALT.find(">") != string::npos) {
        static const std::regex reg("SVSIZE=(\\d+)");
        std::smatch match;
        if (!std::regex_search(INFOSTRUCTTMP.ALT, match, reg)) {
            warning_stream() << "No length information in the fourth column, skip this site -> " << INFOSTRUCTTMP.line << endl;
            return false;
        }

        if (INFOSTRUCTTMP.ALT.find("<INS") != string::npos) {
            refLen = INFOSTRUCTTMP.LEN;
            qryLen = std::stoul(match[1].str());
        } else if (INFOSTRUCTTMP.ALT.find("<DEL") != string::npos) {
            refLen = std::stoul(match[1].str());
            qryLen = INFOSTRUCTTMP.LEN;
        } else if (INFOSTRUCTTMP.ALT.find("<DUP") != string::npos) {
            refLen = std::stoul(match[1].str());
            qryLen = refLen * 2;
        } else {
            warning_stream() << "unknown string -> ref_seq:" << INFOSTRUCTTMP.REF << " qey_seq:" << INFOSTRUCTTMP.ALT << endl;
            refLen = INFOSTRUCTTMP.LEN;

            qryLen = std::stoul(match[1].str());
        }
    } else {
        refLen = INFOSTRUCTTMP.LEN;

        string qrySeq;

        uint32_t maxGtNum = stoi(*max_element(evaluate_gt.begin(), evaluate_gt.end()));
        if (maxGtNum == 0) {  // haploid reference genotype (GT=0)
            return false;
        }
        if (maxGtNum > INFOSTRUCTTMP.ALTVec.size()) {
            throw runtime_error("genotyping and ALT sequence numbers do not match -> " + INFOSTRUCTTMP.line);
        } else {
            qrySeq = INFOSTRUCTTMP.ALTVec[maxGtNum-1];
        }

        qryLen = qrySeq.size();
    }

    if ((INFOSTRUCTTMP.ID.find("Duplication") != string::npos) && (refLen <= 2)) {
        refLen = qryLen;
        qryLen *= 2;
    }

    return true;
}

// Write SNP+indel+ins+del, inversion, duplication and other records to <prefix>.*.vcf.gz; workers measure the records.
void VCFSplit::split_by_type()
{
    VCFReader<SplitResult> VCFOPENCLASS(vcfFileName_, threads_, [this](VCFINFOSTRUCT& INFOSTRUCTTMP, SplitResult& result) {
        measure_record(INFOSTRUCTTMP, result);
    });

    stringstream simpleSVsStream, invSVsStream, dupSVsStream, otherSVsStream;
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10;
    simpleSVsStream.str().reserve(CACHE_SIZE);
    invSVsStream.str().reserve(CACHE_SIZE);
    dupSVsStream.str().reserve(CACHE_SIZE);
    otherSVsStream.str().reserve(CACHE_SIZE);

    SAVE SimpleSave(prefix_ + ".snp.indel.ins.del.vcf.gz", threads_);
    SAVE InvSave(prefix_ + ".inv.vcf.gz", threads_);
    SAVE DupSave(prefix_ + ".dup.vcf.gz", threads_);
    SAVE OtherSave(prefix_ + ".other.vcf.gz", threads_);

    // ------ Classify records by REF and ALT length, flush every 10 MB ------
    ProcessedRecord<SplitResult> record;
    while (VCFOPENCLASS.read(record)) {
        const VCFINFOSTRUCT& INFOSTRUCTTMP = record.info;
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            string headLine = INFOSTRUCTTMP.line + "\n";
            SimpleSave.save(headLine);
            InvSave.save(headLine);
            DupSave.save(headLine);
            OtherSave.save(headLine);
        } else {
            const SplitResult& result = record.get();
            if (!result.keep) {
                continue;
            }
            const uint32_t refLen = result.refLen;
            const uint32_t qryLen = result.qryLen;

            if ((refLen < 50 && qryLen < 50) || (refLen >= 50 && qryLen < 50) || (refLen < 50 && qryLen >= 50)) {
                simpleSVsStream << INFOSTRUCTTMP.line + "\n";
            } else if (refLen >= 50 && qryLen >= 50) {
                if (refLen*2 <= qryLen+10 && refLen*2 >= qryLen-10) {
                    dupSVsStream << INFOSTRUCTTMP.line + "\n";
                } else if (refLen <= qryLen+10 && refLen >= qryLen-10) {
                    invSVsStream << INFOSTRUCTTMP.line + "\n";
                } else {
                    otherSVsStream << INFOSTRUCTTMP.line + "\n";
                }
            } else {
                otherSVsStream << INFOSTRUCTTMP.line + "\n";
            }

            if (simpleSVsStream.tellp() >= CACHE_SIZE ||
                invSVsStream.tellp() >= CACHE_SIZE ||
                dupSVsStream.tellp() >= CACHE_SIZE ||
                otherSVsStream.tellp() >= CACHE_SIZE
            ) {
                string simpleSVs = simpleSVsStream.str();
                string invSVs = invSVsStream.str();
                string dupSVs = dupSVsStream.str();
                string otherSVs = otherSVsStream.str();

                SimpleSave.save(simpleSVs);
                InvSave.save(invSVs);
                DupSave.save(dupSVs);
                OtherSave.save(otherSVs);

                simpleSVsStream.str(string());
                simpleSVsStream.clear();
                invSVsStream.str(string());
                invSVsStream.clear();
                dupSVsStream.str(string());
                dupSVsStream.clear();
                otherSVsStream.str(string());
                otherSVsStream.clear();
            }
        }
    }

    if (simpleSVsStream.tellp() > 0 ||
        invSVsStream.tellp() > 0 ||
        dupSVsStream.tellp() > 0 ||
        otherSVsStream.tellp() > 0
    ) {
        string simpleSVs = simpleSVsStream.str();
        string invSVs = invSVsStream.str();
        string dupSVs = dupSVsStream.str();
        string otherSVs = otherSVsStream.str();

        SimpleSave.save(simpleSVs);
        InvSave.save(invSVs);
        DupSave.save(dupSVs);
        OtherSave.save(otherSVs);

        simpleSVsStream.str(string());
        simpleSVsStream.clear();
        invSVsStream.str(string());
        invSVsStream.clear();
        dupSVsStream.str(string());
        dupSVsStream.clear();
        otherSVsStream.str(string());
        otherSVsStream.clear();
    }
}

// Write each SV to <prefix>.<N>.vcf.gz by the number N of SNPs and indels within -l bp of its breakpoints; workers measure the records.
void VCFSplit::split_by_number()
{
    if (baseVcfFileName_.empty()) {
        baseVcfFileName_ = vcfFileName_;
    }

    stringstream SVs0Stream, SVs1Stream, SVs2Stream, SVs4Stream, SVs6Stream, SVs8Stream, SVsMoreStream;
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10;
    SVs0Stream.str().reserve(CACHE_SIZE);
    SVs1Stream.str().reserve(CACHE_SIZE);
    SVs2Stream.str().reserve(CACHE_SIZE);
    SVs4Stream.str().reserve(CACHE_SIZE);
    SVs6Stream.str().reserve(CACHE_SIZE);
    SVs8Stream.str().reserve(CACHE_SIZE);
    SVsMoreStream.str().reserve(CACHE_SIZE);

    SAVE SVs0Save(prefix_ + ".0.vcf.gz", threads_);
    SAVE SVs1Save(prefix_ + ".1.vcf.gz", threads_);
    SAVE SVs2Save(prefix_ + ".2.vcf.gz", threads_);
    SAVE SVs4Save(prefix_ + ".4.vcf.gz", threads_);
    SAVE SVs6Save(prefix_ + ".6.vcf.gz", threads_);
    SAVE SVs8Save(prefix_ + ".8.vcf.gz", threads_);
    SAVE SVsMoreSave(prefix_ + ".more.vcf.gz", threads_);

    // ------ SNP and indel positions from the base VCF ------
    VCFReader<SplitResult> VCFOPENCLASSBase(baseVcfFileName_, threads_, [this](VCFINFOSTRUCT& INFOSTRUCTTMP, SplitResult& result) {
        measure_record(INFOSTRUCTTMP, result);
    });

    unordered_map<string, snpIndelInfo> snpIndelInfoMap;

    ProcessedRecord<SplitResult> record;
    while (VCFOPENCLASSBase.read(record)) {
        const VCFINFOSTRUCT& INFOSTRUCTBase = record.info;
        if (INFOSTRUCTBase.line.empty()) {
            continue;
        }

        if (INFOSTRUCTBase.line.find("#") == string::npos) {
            const SplitResult& result = record.get();
            if (!result.keep) {
                continue;
            }
            const uint32_t refLen = result.refLen;
            const uint32_t qryLen = result.qryLen;

            uint32_t refEnd;
            refEnd = INFOSTRUCTBase.POS + refLen - 1;

            if ((refLen < 50 && qryLen < 50)) {
                snpIndelInfoMap[INFOSTRUCTBase.CHROM].refStartVec.push_back(INFOSTRUCTBase.POS);
                snpIndelInfoMap[INFOSTRUCTBase.CHROM].refEndVec.push_back(refEnd);
            }
        }
    }

    // ------ SVs of the input VCF ------
    VCFReader<SplitResult> VCFOPENCLASS(vcfFileName_, threads_, [this](VCFINFOSTRUCT& INFOSTRUCTTMP, SplitResult& result) {
        measure_record(INFOSTRUCTTMP, result);
    });

    map<string, map<uint32_t, svInfo> > svInfoMap;

    while (VCFOPENCLASS.read(record)) {
        const VCFINFOSTRUCT& INFOSTRUCTTMP = record.info;
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            string headLine = INFOSTRUCTTMP.line + "\n";
            SVs0Save.save(headLine);
            SVs1Save.save(headLine);
            SVs2Save.save(headLine);
            SVs4Save.save(headLine);
            SVs6Save.save(headLine);
            SVs8Save.save(headLine);
            SVsMoreSave.save(headLine);
        } else {
            const SplitResult& result = record.get();
            if (!result.keep) {
                continue;
            }
            const uint32_t refLen = result.refLen;
            const uint32_t qryLen = result.qryLen;

            uint32_t refEnd;
            refEnd = INFOSTRUCTTMP.POS + refLen - 1;

            if (refLen >= 50 || qryLen >= 50) {
                svInfoMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS].information = INFOSTRUCTTMP.line;
                svInfoMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS].refEnd = refEnd;
            }
        }
    }

    // ------ Count SNPs and indels near both breakpoints of each SV, flush every 10 MB ------
    for (const auto& [chromosome, refStartInfoMap] : svInfoMap) {
        const vector<uint32_t>& snpIndelRefStartVec = snpIndelInfoMap[chromosome].refStartVec;
        const vector<uint32_t>& snpIndelRefEndVec = snpIndelInfoMap[chromosome].refEndVec;

        for (const auto& [refStart, Info] : refStartInfoMap) {
            int64_t startLeftIdxTmp = search_Binary_right(snpIndelRefEndVec, refStart - length_);
            int64_t startRightIdxTmp = search_Binary_left(snpIndelRefEndVec, refStart);
            int64_t endLeftIdxTmp = search_Binary_right(snpIndelRefStartVec, Info.refEnd);
            int64_t endRightIdxTmp = search_Binary_left(snpIndelRefStartVec, Info.refEnd + length_);

            int64_t leftSnpIndelNum = startRightIdxTmp - startLeftIdxTmp;
            if (leftSnpIndelNum < 0) {
                leftSnpIndelNum = 0;
            } else {
                leftSnpIndelNum++;
            }

            int64_t rightSnpIndelNum = endRightIdxTmp - endLeftIdxTmp;
            if (rightSnpIndelNum < 0) {
                rightSnpIndelNum = 0;
            } else {
                rightSnpIndelNum++;
            }

            int64_t snpIndelNum = leftSnpIndelNum + rightSnpIndelNum;
            if (snpIndelNum == 0) {
                SVs0Stream << Info.information + "\n";
            } else if (snpIndelNum == 1) {
                SVs1Stream << Info.information + "\n";
            } else if (snpIndelNum == 2) {
                SVs2Stream << Info.information + "\n";
            } else if (snpIndelNum == 4) {
                SVs4Stream << Info.information + "\n";
            } else if (snpIndelNum == 6) {
                SVs6Stream << Info.information + "\n";
            } else if (snpIndelNum == 8) {
                SVs8Stream << Info.information + "\n";
            } else {
                SVsMoreStream << Info.information + "\n";
            }

            if (SVs0Stream.tellp() >= CACHE_SIZE ||
                SVs1Stream.tellp() >= CACHE_SIZE ||
                SVs2Stream.tellp() >= CACHE_SIZE ||
                SVs4Stream.tellp() >= CACHE_SIZE ||
                SVs6Stream.tellp() >= CACHE_SIZE ||
                SVs8Stream.tellp() >= CACHE_SIZE ||
                SVsMoreStream.tellp() >= CACHE_SIZE
            ) {
                string SVs0 = SVs0Stream.str();
                string SVs1 = SVs1Stream.str();
                string SVs2 = SVs2Stream.str();
                string SVs4 = SVs4Stream.str();
                string SVs6 = SVs6Stream.str();
                string SVs8 = SVs8Stream.str();
                string SVsMore = SVsMoreStream.str();

                SVs0Save.save(SVs0);
                SVs1Save.save(SVs1);
                SVs2Save.save(SVs2);
                SVs4Save.save(SVs4);
                SVs6Save.save(SVs6);
                SVs8Save.save(SVs8);
                SVsMoreSave.save(SVsMore);

                SVs0Stream.str(string());
                SVs0Stream.clear();
                SVs1Stream.str(string());
                SVs1Stream.clear();
                SVs2Stream.str(string());
                SVs2Stream.clear();
                SVs4Stream.str(string());
                SVs4Stream.clear();
                SVs6Stream.str(string());
                SVs6Stream.clear();
                SVs8Stream.str(string());
                SVs8Stream.clear();
                SVsMoreStream.str(string());
                SVsMoreStream.clear();
            }
        }
    }

    if (SVs0Stream.tellp() > 0 ||
        SVs1Stream.tellp() > 0 ||
        SVs2Stream.tellp() > 0 ||
        SVs4Stream.tellp() > 0 ||
        SVs6Stream.tellp() > 0 ||
        SVs8Stream.tellp() > 0 ||
        SVsMoreStream.tellp() > 0
    ) {
        string SVs0 = SVs0Stream.str();
        string SVs1 = SVs1Stream.str();
        string SVs2 = SVs2Stream.str();
        string SVs4 = SVs4Stream.str();
        string SVs6 = SVs6Stream.str();
        string SVs8 = SVs8Stream.str();
        string SVsMore = SVsMoreStream.str();

        SVs0Save.save(SVs0);
        SVs1Save.save(SVs1);
        SVs2Save.save(SVs2);
        SVs4Save.save(SVs4);
        SVs6Save.save(SVs6);
        SVs8Save.save(SVs8);
        SVsMoreSave.save(SVsMore);

        SVs0Stream.str(string());
        SVs0Stream.clear();
        SVs1Stream.str(string());
        SVs1Stream.clear();
        SVs2Stream.str(string());
        SVs2Stream.clear();
        SVs4Stream.str(string());
        SVs4Stream.clear();
        SVs6Stream.str(string());
        SVs6Stream.clear();
        SVs8Stream.str(string());
        SVs8Stream.clear();
        SVsMoreStream.str(string());
        SVsMoreStream.clear();
    }
}
