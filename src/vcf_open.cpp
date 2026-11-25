#include "../include/vcf_open.hpp"

#include <functional>
#include <regex>
#include <stdexcept>

#include "../include/logger.hpp"
#include "../include/vcf_sort.hpp"

using namespace std;

namespace {

// Split a GT value into its allele tokens at every '/' and '|'.
vector<string> split_alleles(const string & gtTxt) {
    vector<string> alleles;
    string allele;
    for (const char c : gtTxt) {
        if (c == '/' || c == '|') {
            alleles.push_back(allele);
            allele.clear();
        } else {
            allele += c;
        }
    }
    alleles.push_back(allele);
    return alleles;
}

// Return true if every allele of a GT value is '.' (any ploidy).
bool is_missing_gt(const string & gtTxt) {
    for (const auto& allele : split_alleles(gtTxt)) {
        if (allele != ".") {
            return false;
        }
    }
    return true;
}

// Return the end of the last complete record in text (just after a '\n' that closes a non-blank line), looking only at '\n's from `from` on; 0 if there is none.
size_t find_record_end(const string & text, size_t from) {
    for (size_t pos = text.size(); pos > from; ) {
        const size_t newline = text.rfind('\n', pos - 1);
        if (newline == string::npos || newline < from) break;
        if (newline > 0 && text[newline - 1] != '\n') return newline + 1;
        pos = newline;
    }
    return 0;
}

}

// Open a plain, gzip or BGZF VCF file; throw if it cannot be opened.
VCFFile::VCFFile(
    const string & vcfFileName,
    uint32_t threads
) : vcfFileName_(vcfFileName), threads_(threads), file_(vcfFileName, threads), eof_(false) {}

// Read the next ~256 KB of whole records (the rest at EOF); false at the end; throw on damaged input or a missing BGZF EOF marker.
bool VCFFile::read_chunk(
    string & chunk
) {
    static constexpr size_t chunkSize = 256 << 10;

    size_t end = 0;
    while (end == 0 && !eof_) {
        const size_t from = pending_.size();
        pending_.resize(from + chunkSize);
        const ssize_t bytes = file_.read(&pending_[from], chunkSize);
        pending_.resize(from + bytes);
        eof_ = bytes == 0;
        end = eof_ ? pending_.size() : find_record_end(pending_, from);
    }

    chunk.swap(pending_);
    pending_.assign(chunk, end, string::npos);
    chunk.resize(end);
    return !chunk.empty();
}

// Split a chunk into record lines (views; '\n' dropped, '\r' kept, a blank line joined to the next, so a final blank line is empty).
vector<string_view> VCFFile::split_lines(
    const string & chunk
) {
    vector<string_view> lines;
    size_t pos = 0;
    while (pos < chunk.size()) {
        if (chunk[pos] == '\n' && ++pos == chunk.size()) {
            lines.emplace_back();
            break;
        }
        size_t end = chunk.find('\n', pos);
        if (end == string::npos) end = chunk.size();
        lines.emplace_back(chunk.data() + pos, end - pos);
        pos = end + 1;
    }
    return lines;
}

// Fill INFOSTRUCTTMP from one record line, reusing its strings: lines containing '#' keep only line and lineVec; other lines need at least 9 columns.
void VCFFile::parse_record(
    string_view line,
    VCFINFOSTRUCT & INFOSTRUCTTMP
) const {
    if (line.find('#') != string_view::npos) {
        INFOSTRUCTTMP.clear();
        INFOSTRUCTTMP.line.assign(line.data(), line.size());
        split(INFOSTRUCTTMP.line, "\t", INFOSTRUCTTMP.lineVec);
        return;
    }

    INFOSTRUCTTMP.line.assign(line.data(), line.size());
    split(INFOSTRUCTTMP.line, "\t", INFOSTRUCTTMP.lineVec);

    if (INFOSTRUCTTMP.lineVec.size() < 9) {
        throw runtime_error("'" + vcfFileName_ + "' -> The number of vcf columns is less than 9.\n" + INFOSTRUCTTMP.line);
    }

    INFOSTRUCTTMP.CHROM = INFOSTRUCTTMP.lineVec[0];
    INFOSTRUCTTMP.POS = stoul(INFOSTRUCTTMP.lineVec[1]);
    INFOSTRUCTTMP.ID = INFOSTRUCTTMP.lineVec[2];
    INFOSTRUCTTMP.REF = INFOSTRUCTTMP.lineVec[3];
    INFOSTRUCTTMP.ALT = INFOSTRUCTTMP.lineVec[4];
    split(INFOSTRUCTTMP.ALT, ",", INFOSTRUCTTMP.ALTVec);

    if (isdigit(INFOSTRUCTTMP.lineVec[5][0])) {
        INFOSTRUCTTMP.QUAL = stod(INFOSTRUCTTMP.lineVec[5]);
    } else {
        INFOSTRUCTTMP.QUAL = 0.0;
    }

    INFOSTRUCTTMP.FILTER = INFOSTRUCTTMP.lineVec[6];
    INFOSTRUCTTMP.INFO = INFOSTRUCTTMP.lineVec[7];
    INFOSTRUCTTMP.FORMAT = INFOSTRUCTTMP.lineVec[8];

    INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.REF.size();
    INFOSTRUCTTMP.END = INFOSTRUCTTMP.POS + INFOSTRUCTTMP.LEN - 1;
}

// Classify a record as SNP, InDel, Inversion, Duplication, Deletion, Insertion or Other from its allele lengths.
string VCFFile::get_TYPE(
    const uint32_t & refLen,
    const vector<string> & ALTVec
) {
    uint32_t qryLen = 0;
    for (const auto& it1 : ALTVec) {
        uint32_t qryLenTmp = it1.size();
        qryLen = max(qryLen, qryLenTmp);
    }

    int32_t svLen = qryLen - refLen;
    double lengthRatio = qryLen / float(refLen);

    if (svLen == 0 && refLen == 1 && qryLen == 1) {
        return "SNP";
    } else if (svLen <= 49 && svLen >= -49 && refLen <= 49 && qryLen <= 49) {
        return "InDel";
    } else if (svLen >= -2 && svLen <= 2 && refLen > 49 && qryLen > 49) {
        return "Inversion";
    } else if (lengthRatio >= 1.8 && lengthRatio <= 2.2 && refLen > 49 && qryLen > 49) {
        return "Duplication";
    } else if (svLen < 0) {
        return "Deletion";
    } else if (svLen > 0) {
        return "Insertion";
    } else {
        return "Other";
    }
}

// Return the genotype alleles of every sample ('.' read as '0') and the number of missing samples (all alleles '.', or no GT field).
tuple<unordered_map<int, vector<string> >, uint32_t> VCFFile::get_gt(
    const vector<string> & lineTmpVec
) {
    unordered_map<int, vector<string> > GTVecMap;
    uint32_t misSampleNum = 0;

    vector<string> formatVec = split(lineTmpVec[8], ":");

    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));

    if (gtIndex == formatVec.size()) {
        warning_stream() << "no [GT] information in FORMAT -> " << lineTmpVec[0] << ":" << lineTmpVec[1] << endl;
        for (size_t i = 9; i < lineTmpVec.size(); i++) {
            GTVecMap[i - 9] = {"0", "0"};
            misSampleNum++;
        }
    } else {
        string GTString;

        for (size_t i = 9; i < lineTmpVec.size(); i++) {
            GTString = sample_gt(lineTmpVec, i, gtIndex);

            if (is_missing_gt(GTString)) {
                misSampleNum++;
            }

            vector<string> GTVec = gt_split(GTString);
            GTVecMap[i - 9] = GTVec;
        }
    }

    return tuple(GTVecMap, misSampleNum);
}

// Split a genotype into alleles, replacing '.' with '0'; throw if it is not a number or separated by '/' or '|'.
vector<string> VCFFile::gt_split(
    string & gtTxt
) {
    vector<string> gtVecTmp;

    std::replace(gtTxt.begin(), gtTxt.end(), '.', '0');

    if (gtTxt.find("/") != string::npos) {
        gtVecTmp = split(gtTxt, "/");
    } else if (gtTxt.find("|") != string::npos) {
        gtVecTmp = split(gtTxt, "|");
    } else {
        try {
            warning_stream() << "sample has only one genotype, attempting to correct to diploid -> " << gtTxt << endl;
            stoul(gtTxt);
            gtVecTmp.push_back(gtTxt);
        } catch (const std::logic_error&) {  // invalid_argument or out_of_range
            throw runtime_error("GT is not separated by '/' or '|' -> " + gtTxt);
        }
    }

    return gtVecTmp;
}

// Return the GT value (field gtIndex) of a sample column; throw if the column is missing or has fewer fields.
string VCFFile::sample_gt(
    const vector<string> & lineVec,
    size_t column,
    uint32_t gtIndex
) {
    if (column >= lineVec.size()) {
        throw runtime_error("The sample column is missing -> " + lineVec[0] + " " + lineVec[1]);
    }
    vector<string> fields = split(lineVec[column], ":");
    if (gtIndex >= fields.size()) {
        throw runtime_error("The sample column has no genotype (GT) field -> " + lineVec[0] + " " + lineVec[1] + " " + lineVec[column]);
    }
    return std::move(fields[gtIndex]);
}

// Return the minor allele frequency among called alleles (second-largest allele count / called alleles, 0 with one allele) and the missing rate (missing samples / samples).
tuple<double, double> VCFFile::calculate(
    const vector<string> & lineVec,
    uint32_t misSampleNum
) {
    double MAFTMP = 0;
    double MISSRATETMP;

    uint32_t sampleNum = lineVec.size() - 9;

    vector<string> formatVec = split(lineVec[8], ":");
    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));

    map<string, uint32_t> alleleCountMap;
    uint32_t calledAlleleNum = 0;
    if (gtIndex < formatVec.size()) {
        for (size_t i = 9; i < lineVec.size(); i++) {
            for (const auto& allele : split_alleles(sample_gt(lineVec, i, gtIndex))) {
                if (allele == ".") {
                    continue;
                }
                alleleCountMap[allele]++;
                calledAlleleNum++;
            }
        }
    }

    vector<uint32_t> alleleCountVec;
    for (const auto& it1 : alleleCountMap) {
        alleleCountVec.push_back(it1.second);
    }
    sort(alleleCountVec.begin(), alleleCountVec.end(), greater<uint32_t>());

    if (alleleCountVec.size() > 1) {
        MAFTMP = alleleCountVec[1]/(double)calledAlleleNum;
    }

    MISSRATETMP = misSampleNum/(double)sampleNum;

    return make_tuple(MAFTMP, MISSRATETMP);
}

// Throw if the records of a VCF file are not sorted by position within each chromosome (the order rule of SortedIndexPass).
void check_vcf_sort(const string & inputVcf)
{
    SortedIndexPass pass(inputVcf);
    try {
        VCFINFOSTRUCT INFOSTRUCTTMP;
        VCFOPEN VCFOPENCLASS(inputVcf);

        while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
            if (INFOSTRUCTTMP.line.empty() || INFOSTRUCTTMP.line.find("#") != string::npos) {
                continue;
            }
            pass.check_order(INFOSTRUCTTMP);
        }
    } catch (...) {
        pass.fail<bool>(current_exception()).get();
    }
    pass.finish(true).get();
}

// Return the REF length and the ALT lengths of a record: every ALT allele ("all") or the alleles in gtVec ("hap").
tuple<uint32_t, vector<uint32_t> > get_hap_len(
    const string & svType,
    const string & refSeq,
    const string & qrySeqs,
    const vector<int> & gtVec,
    const string & lenType
)
{
    if (lenType != "hap" && lenType != "all") {
        throw runtime_error("lenType -> " + lenType);
    }

    uint32_t refLen = refSeq.size();
    vector<string> qrySeqVec = split(qrySeqs, ",");

    if (lenType == "hap") {
        int maxGtNum = *max_element(gtVec.begin(), gtVec.end());
        if (static_cast<uint32_t>(maxGtNum) > qrySeqVec.size()) {
            throw runtime_error("number of genotyping and ALT sequences do not match -> " + qrySeqs);
        }
    }

    // ------ Length of every allele (index 0 = REF) ------
    vector<uint32_t> seqLenVec;
    seqLenVec.push_back(refLen);

    for (size_t i = 0; i < qrySeqVec.size(); i++) {
        const string& qrySeq = qrySeqVec[i];

        uint32_t refLenTmp = refLen;
        uint32_t qryLenTmp = qrySeq.size();

        if (qrySeq.find(">") != string::npos) {
            static const std::regex reg("SVSIZE=(\\d+)");
            std::smatch match;
            if (std::regex_search(qrySeq, match, reg)) {
                if (qrySeq.find("<INS") != string::npos && qrySeq.find("<DEL") == string::npos) {
                    refLenTmp = refLen;
                    qryLenTmp = std::stoul(match[1].str());;
                }
                else if (qrySeq.find("<INS") == string::npos && qrySeq.find("<DEL") != string::npos) {
                    refLenTmp = std::stoul(match[1].str());;
                    qryLenTmp = refLen;
                }
                else if (qrySeq.find("<DUP") != string::npos) {
                    refLenTmp = std::stoul(match[1].str());;
                    qryLenTmp = refLenTmp * 2;
                }
            } else {
                warning_stream() << "No length information in the fourth column -> " << qrySeq << endl;
                refLenTmp = refLen;
                qryLenTmp = qrySeq.size();
            }

            seqLenVec[0] = refLenTmp;
        }

        if (svType.find("Duplication") != string::npos && refLenTmp <= 2) {
            refLenTmp = qryLenTmp;
            qryLenTmp *= 2;
        }

        if (seqLenVec.size() == (i + 1)) {
            seqLenVec.push_back(qryLenTmp);
        } else {
            throw runtime_error("Incorrect index for haplotype length -> " + qrySeqs);
        }
    }

    // ------ Select the requested alleles ------
    vector<uint32_t> qryLenVec;
    if (lenType == "hap") {
        for (auto gtTmp : gtVec) {
            qryLenVec.push_back(seqLenVec[gtTmp]);
        }
        if (gtVec.size() == 1) {
            qryLenVec.push_back(qryLenVec[0]);  // a haploid genotype counts its allele twice
        }
    }
    else {
        for (size_t i = 1; i < seqLenVec.size(); i++) {
            qryLenVec.push_back(seqLenVec[i]);
        }
    }

    return make_tuple(seqLenVec[0], qryLenVec);
}
