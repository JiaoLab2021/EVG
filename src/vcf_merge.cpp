#include "../include/vcf_merge.hpp"

#include <deque>
#include <future>
#include <stdexcept>
#include <utility>

#include "../include/logger.hpp"
#include "../include/task_pool.hpp"
#include "../include/vcf_sort.hpp"

using namespace std;

// Store the input files, sample name, mode, output name and thread count.
VCFMerge::VCFMerge(
    const string& trueVcf,
    const string& ParagraphVcf,
    const string& GraphTyper2Vcf,
    const string& BayesTyperVcf,
    const string& MAPVcf,
    const string& GiraffeVcf,
    const string& GraphAlignerVcf,
    const string& PanGenieVcf,
    const string& sampleName,
    const string& mode,
    const string& outputFileName,
    uint32_t threads
) : trueVcf_(trueVcf), ParagraphVcf_(ParagraphVcf), GraphTyper2Vcf_(GraphTyper2Vcf),
    BayesTyperVcf_(BayesTyperVcf), MAPVcf_(MAPVcf), GiraffeVcf_(GiraffeVcf),
    GraphAlignerVcf_(GraphAlignerVcf), PanGenieVcf_(PanGenieVcf), sampleName_(sampleName),
    threads_(threads), mMode(mode), mOutputFileName(outputFileName) {}

// Index the population VCF (output header, samples, positions, allele lengths, genotypes) and check its order in the same pass.
Outcome<baseVcfStruct> VCFMerge::build_basefile_index() const {
    baseVcfStruct index;

    log_stream() << "Building index: '" << trueVcf_ << "'" << endl;

    SortedIndexPass pass(trueVcf_);
    try {
        VCFINFOSTRUCT INFOSTRUCTTMP;
        VCFOPEN VCFOPENCLASS(trueVcf_, threads_);

        while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
            if (INFOSTRUCTTMP.line.find("#") != string::npos) {
                pass.index([&] { index_base_header(INFOSTRUCTTMP, index); });
            } else {
                pass.check_order(INFOSTRUCTTMP);
                pass.index([&] { index_base_record(INFOSTRUCTTMP, index); });
            }
        }
    } catch (...) {
        return pass.fail<baseVcfStruct>(current_exception());
    }

    return pass.finish(std::move(index));
}

// Add a header line of the population VCF to the output header; #CHROM also gives the sample names and the column count.
void VCFMerge::index_base_header(const VCFINFOSTRUCT & INFOSTRUCTTMP, baseVcfStruct & index) const {
    if (INFOSTRUCTTMP.line.find("#CHROM") != string::npos) {
        index.headInfo += "##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Software\">\n";
        index.headInfo += "##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Total filtered read depth used for genotyping.\">\n";
        index.headInfo += "##FORMAT=<ID=ADP,Number=1,Type=Float,Description=\"Average depth\">\n";
        index.headInfo += "##FORMAT=<ID=NDP,Number=1,Type=Float,Description=\"Depth after z-score normalization\">\n";

        vector<string> informationsVecTmp{&INFOSTRUCTTMP.lineVec[0], &INFOSTRUCTTMP.lineVec[0]+9};
        index.headInfo += join(informationsVecTmp, "\t") + "\t" + sampleName_ + "\n";

        vector<string> sampleNameVecTmp{&INFOSTRUCTTMP.lineVec[9], &INFOSTRUCTTMP.lineVec[0] + INFOSTRUCTTMP.lineVec.size()};
        index.sampleNameVec = sampleNameVecTmp;

        index.colNum = INFOSTRUCTTMP.lineVec.size();
    } else {
        index.headInfo += INFOSTRUCTTMP.line + "\n";
    }
}

// Index one population record: position, site columns with genotypes, REF length and every ALT length.
void VCFMerge::index_base_record(const VCFINFOSTRUCT & INFOSTRUCTTMP, baseVcfStruct & index) const {
    index.startMap[INFOSTRUCTTMP.CHROM].push_back(INFOSTRUCTTMP.POS);

    vector<string> gtVec = this->get_gt_all(INFOSTRUCTTMP.lineVec);

    vector<string> informationsVecTmp{&INFOSTRUCTTMP.lineVec[0], &INFOSTRUCTTMP.lineVec[0]+8};
    informationsVecTmp.push_back("GT:SO:DP:ADP:NDP");

    index.BaseInfoTupMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS] = make_tuple(join(informationsVecTmp, "\t"), gtVec);

    uint32_t refLen;
    vector<uint32_t> qryLenVec;
    tie(refLen, qryLenVec) = get_hap_len(
        INFOSTRUCTTMP.ID,
        INFOSTRUCTTMP.REF,
        INFOSTRUCTTMP.ALT,
        {0, 0},
        "all"
    );

    index.refLenMap[INFOSTRUCTTMP.CHROM].push_back(refLen);
    index.qryLenVecVecMap[INFOSTRUCTTMP.CHROM].push_back(qryLenVec);
}

// Return the genotypes of all samples with '|' separators; missing or unknown genotypes become 0|0.
vector<string> VCFMerge::get_gt_all(
    const vector<string> & lineVec
) const {
    vector<string> gtVec;

    vector<string> formatVec = split(lineVec[8], ":");

    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));

    if (gtIndex == formatVec.size()) {
        warning_stream() << "No [GT] information in FORMAT, replace with 0|0 -> " << lineVec[0] << ":" << lineVec[1] << endl;
        for (size_t i = 9; i < lineVec.size(); i++) {
            gtVec.push_back("0|0");
        }
    } else {
        for (size_t i = 9; i < lineVec.size(); i++) {
            string gt = VCFOPEN::sample_gt(lineVec, i, gtIndex);

            string splitStr;
            if (gt.find("/") != string::npos) {
                splitStr = "/";
            } else if (gt.find("|") != string::npos) {
                splitStr = "|";
            } else {
                warning_stream() << "Genotype does not contain '\\' or '|', replace with 0|0 -> " << lineVec[0] << ":" << lineVec[1] << endl;
                gtVec.push_back("0|0");
                continue;
            }

            if (gt == "." || gt == ".|." || gt == "./." || gt == "0" || gt == "0|0" || gt == "0/0") {
                gtVec.push_back("0|0");
                continue;
            }

            replace(gt.begin(), gt.end(), '/', '|');
            gtVec.push_back(gt);
        }
    }

    return gtVec;
}

// Index the population VCF here while up to threads - 1 tool VCFs are indexed as tasks (-t 1: one by one after it); merge each tool in order, then release it.
void VCFMerge::run_index_merge()
{
    vector<string> softwareVec = {ParagraphVcf_, MAPVcf_, GraphTyper2Vcf_, BayesTyperVcf_, GiraffeVcf_, GraphAlignerVcf_, PanGenieVcf_};
    vector<string> softwareNameVec = {"Paragraph", "VG-MAP", "GraphTyper2", "BayesTyper", "VG-Giraffe", "GraphAligner", "PanGenie"};

    vector<size_t> toolIdx;
    for (size_t i = 0; i < softwareVec.size(); i++) {
        if (!softwareVec[i].empty()) {
            toolIdx.push_back(i);
        }
    }

    const size_t maxTasks = threads_ > 1 ? threads_ - 1 : 1;
    deque<pair<size_t, future<Outcome<softwareVcfStruct> > > > tasks;  // tool index tasks in tool order (each slot exists before its task)
    size_t nextTool = 0;

    // Start tool index tasks up to the budget.
    auto start_tools = [&] {
        while (nextTool < toolIdx.size() && tasks.size() < maxTasks) {
            const size_t i = toolIdx[nextTool++];
            tasks.emplace_back(i, future<Outcome<softwareVcfStruct> >());
            tasks.back().second = start_task(threads_, [this, file = softwareVec[i], software = softwareNameVec[i]] {
                return build_softwarefile_index(file, software);
            });
        }
    };

    try {
        if (threads_ > 1) start_tools();

        mergeVcfStruct_ = std::move(build_basefile_index().get());

        // ------ Commit in tool order: messages, then the index or the error, then the merge ------
        while (true) {
            start_tools();
            if (tasks.empty()) break;

            const size_t i = tasks.front().first;
            softwareVcfStruct softvcfStructure = std::move(tasks.front().second.get().get());
            tasks.pop_front();

            mergeVcfStruct_.softwateSampleIdx[softwareNameVec[i]] = softvcfStructure.sampleIdx;
            mergeVcfStruct_.depthMap[softwareNameVec[i]] = softvcfStructure.depthStats;

            vcf_merge(softvcfStructure);

            softwareNum_++;

            if (softwareNameVec[i] == "Paragraph") {
                paragraphInputBool_ = true;
            }
        }
    } catch (...) {
        for (auto& task : tasks) {
            if (task.second.valid()) task.second.wait();
        }
        throw;
    }
}

// Index the PASS, non-reference calls, sample column and depth statistics of one tool VCF, checking its order (a pool task: one reading thread).
Outcome<softwareVcfStruct> VCFMerge::build_softwarefile_index(
    const string & vcfFileName,
    const string & software
) const
{
    softwareVcfStruct softvcfStructure;

    softvcfStructure.software = software;

    SortedIndexPass pass(vcfFileName);
    try {
        log_stream() << "Building index: " << software << endl;

        VCFINFOSTRUCT INFOSTRUCTTMP;
        VCFOPEN VCFOPENCLASS(vcfFileName);

        while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
            if (INFOSTRUCTTMP.line.find("#") != string::npos) {
                pass.index([&] { index_software_header(INFOSTRUCTTMP, softvcfStructure); });
            } else {
                pass.check_order(INFOSTRUCTTMP);
                pass.index([&] { index_software_record(INFOSTRUCTTMP, softvcfStructure); });
            }
        }
    } catch (...) {
        return pass.fail<softwareVcfStruct>(current_exception());
    }

    // ------ Depth statistics (not after an index error) ------
    pass.index([&] { softvcfStructure.depthStats = cal_var_sd(softvcfStructure.depthVec); });

    return pass.finish(std::move(softvcfStructure));
}

// Find the sample's column on the #CHROM line of a tool VCF; throw if the sample is not there.
void VCFMerge::index_software_header(const VCFINFOSTRUCT & INFOSTRUCTTMP, softwareVcfStruct & softvcfStructure) const {
    if (INFOSTRUCTTMP.line.find("#CHROM") != string::npos) {
        uint32_t sampleIdx = distance(INFOSTRUCTTMP.lineVec.begin(), find(INFOSTRUCTTMP.lineVec.begin()+9, INFOSTRUCTTMP.lineVec.end(), sampleName_));

        if (sampleIdx == INFOSTRUCTTMP.lineVec.size()) {
            throw runtime_error("'" + sampleName_ + "' is not in the output file of '"
                + softvcfStructure.software + "' -> " + INFOSTRUCTTMP.line);
        }

        softvcfStructure.sampleIdx = sampleIdx;
    }
}

// Index one tool call: its depth always counts for the statistics; PASS, non-reference calls also get position, genotype and lengths.
void VCFMerge::index_software_record(const VCFINFOSTRUCT & INFOSTRUCTTMP, softwareVcfStruct & softvcfStructure) const {
    const uint32_t sampleIdx = softvcfStructure.sampleIdx;
    if (sampleIdx == 0) {
        throw runtime_error("'" + sampleName_ + "' is not in the output file of '" + softvcfStructure.software + "'");
    }
    if (sampleIdx >= INFOSTRUCTTMP.lineVec.size()) {
        throw runtime_error("The sample column is missing -> " + INFOSTRUCTTMP.CHROM + " " + to_string(INFOSTRUCTTMP.POS));
    }

    vector<float> depthVec = get_depth(
        INFOSTRUCTTMP.lineVec,
        sampleIdx
    );
    float depSum = accumulate(depthVec.begin(), depthVec.end(), 0.0);
    softvcfStructure.depthVec.push_back(depSum);

    vector<int> gtVec = get_gt_sample(
        INFOSTRUCTTMP.lineVec,
        sampleIdx
    );
    string gt = join(gtVec, "/");

    if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {
        warning_stream() << "The genotype of the variant is empty, skip this site -> "
            << INFOSTRUCTTMP.CHROM << " "
            << INFOSTRUCTTMP.POS << endl;
        return;
    }

    if (INFOSTRUCTTMP.FILTER != "PASS") {
        return;
    }

    int refLen;
    vector<uint32_t> qryLenVec;
    tie(refLen, qryLenVec) = get_hap_len(
        INFOSTRUCTTMP.ID,
        INFOSTRUCTTMP.REF,
        INFOSTRUCTTMP.ALT,
        gtVec,
        "hap"
    );

    softvcfStructure.startMap[INFOSTRUCTTMP.CHROM].push_back(INFOSTRUCTTMP.POS);
    softvcfStructure.gtVecMap[INFOSTRUCTTMP.CHROM].push_back(gtVec);
    softvcfStructure.depthVecMap[INFOSTRUCTTMP.CHROM].push_back(depthVec);

    softvcfStructure.refLenMap[INFOSTRUCTTMP.CHROM].push_back(refLen);
    softvcfStructure.qryLenVecVecMap[INFOSTRUCTTMP.CHROM].push_back(qryLenVec);
}

// Return the per-allele depth of the sample from FORMAT DP, MAC (BayesTyper) or KC (PanGenie).
vector<float> VCFMerge::get_depth(
    const vector<string> & informationsVec,
    const uint32_t & sampleIdx
) const {
    vector<float> depthVec;

    vector<string> formatVec = split(informationsVec[8], ":");

    uint32_t depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "DP"));

    if (depthIndex == formatVec.size()) {
        depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "MAC"));
    }
    if (depthIndex == formatVec.size()) {
        depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "KC"));
    }

    if (depthIndex == formatVec.size()) {
        warning_stream() << "no [DP,MAC,KC] information in FORMAT -> " << informationsVec[0] << ":" << informationsVec[1] << endl;
        depthVec = {0, 0};
    }
    else {
        vector<string> splitResult = split(informationsVec[sampleIdx], ":");

        if (depthIndex < splitResult.size()) {
            auto depthSplit = split(splitResult[depthIndex], ",");
            for (auto it : depthSplit) {
                if (it == "-1") {
                    it = "0";
                }
                try {
                    depthVec.push_back(stof(it));
                } catch (const std::exception& e) {
                    warning_stream() << "stof(it) conversion failed. Replaced with {0,0} -> " << it << endl;
                    depthVec = {0, 0};
                }
            }
        }
        else {
            warning_stream() << "no [DP,MAC,KC] information in column " << sampleIdx << " -> " << informationsVec[sampleIdx] << endl;
            depthVec = {0, 0};
        }
    }

    return depthVec;
}

// Return the genotype of one sample (0 = last column) as integers; missing or unparsable genotypes become {0, 0}.
vector<int> VCFMerge::get_gt_sample(
    const vector<string> & informationsVec,
    uint32_t sampleIdx
) const {
    vector<int> gtVec;

    vector<string> formatVec = split(informationsVec[8], ":");

    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));

    if (gtIndex == formatVec.size()) {
        warning_stream() << "No [GT] information in FORMAT, replace with 0|0 -> " << informationsVec[0] << ":" << informationsVec[1] << endl;
        gtVec = {0, 0};
    } else {
        string gt;

        if (sampleIdx == 0) {
            gt = VCFOPEN::sample_gt(informationsVec, informationsVec.size() - 1, gtIndex);
        }
        else {
            gt = VCFOPEN::sample_gt(informationsVec, sampleIdx, gtIndex);
        }

        string splitStr;
        if (gt.find("/") != string::npos) {
            splitStr = "/";
        }  else if (gt.find("|") != string::npos) {
            splitStr = "|";
        }  else {
            warning_stream() << "Genotype does not contain '\\' or '|', replace with 0|0 -> " << informationsVec[0] << ":" << informationsVec[1] << endl;
            gtVec = {0, 0};
            return gtVec;
        }

        auto splitResult = split(gt, splitStr);
        for (auto it : splitResult) {
            if (it == ".") {
                gtVec = {0, 0};
                return gtVec;
            }
            try {
                gtVec.push_back(stoi(it));
            } catch (const std::exception& e) {
                warning_stream() << "stoi(it) conversion failed. Replaced with {0,0} -> " << it << endl;
                gtVec = {0, 0};
                return gtVec;
            }
        }
    }

    return gtVec;
}

// Record the depth and genotype of a matched tool allele at a population site; return -1 if only one side is reference.
tuple<int, int> VCFMerge::recall_push(
    const string chromosome,
    const uint32_t trueRefStart,
    const string software,
    const vector<uint32_t> qryLenVec,
    const vector<int> gtVec,
    const int gt,
    const float depth,
    uint32_t j,
    uint32_t k,
    uint32_t l,
    int64_t & indexLeft,
    int64_t & indexRight
)
{
    if (mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome].find(trueRefStart) == mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome].end()) {
        mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome][trueRefStart];
    }
    if (mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome][trueRefStart].find(software) == mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome][trueRefStart].end()) {
        vector<float> depthVecTmp(qryLenVec.size());
        vector<int> gtVecTmp(qryLenVec.size());
        mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software] = make_tuple(depthVecTmp, gtVecTmp);
    }

    if ((gt != 0 && l == 0) || (gt == 0 && l != 0)) {
        return make_tuple(-1, j);
    } else {
        get<0>(mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = depth;
        get<1>(mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = l;

        indexLeft = k;
        indexRight = k;

        if (j == 0 && gt == gtVec[1]) {
            j++;
            get<0>(mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = depth;
            get<1>(mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = l;
        }
    }

    return make_tuple(0, j);
}

// Match every allele of one tool to population sites by position and length tolerance.
void VCFMerge::vcf_merge(
    softwareVcfStruct & softvcfStructure
) {
    vector<string> chromosomeVec;
    string software = softvcfStructure.software;

    log_stream() << "Merging: " << software << ".\n";

    for (auto iter = mergeVcfStruct_.startMap.begin(); iter != mergeVcfStruct_.startMap.end(); iter++) {
        chromosomeVec.push_back(iter->first);
    }

    // ------ Candidate population sites within 10 bp (both alleles < 50 bp) or 200 bp ------
    for (const auto& [chromosome, startVec] : softvcfStructure.startMap) {
        const vector<uint32_t>& refLenVec = softvcfStructure.refLenMap[chromosome];
        const vector<vector<uint32_t> >& qryLenVecVec = softvcfStructure.qryLenVecVecMap[chromosome];
        const vector<vector<int> >& gtVecVec = softvcfStructure.gtVecMap[chromosome];
        const vector<vector<float> >& depthVecVec = softvcfStructure.depthVecMap[chromosome];

        if (mergeVcfStruct_.recallSoftwareGtDepVecMap.find(chromosome) == mergeVcfStruct_.recallSoftwareGtDepVecMap.end()) {
            mergeVcfStruct_.recallSoftwareGtDepVecMap[chromosome];
        }

        if (find(chromosomeVec.begin(), chromosomeVec.end(), chromosome) == chromosomeVec.end()) {
            throw runtime_error("chromosome is not present in '" + trueVcf_ + "' -> " + chromosome);
        } else {
            const vector<uint32_t>& mergeRefStartVec = mergeVcfStruct_.startMap[chromosome];
            const vector<uint32_t>& mergeRefLenVec = mergeVcfStruct_.refLenMap[chromosome];
            const vector<vector<uint32_t> >& mergeQryLenVecVec = mergeVcfStruct_.qryLenVecVecMap[chromosome];

            for (size_t i = 0; i < startVec.size(); i++) {
                uint32_t refStart = startVec[i];
                uint32_t refLen = refLenVec[i];
                vector<uint32_t> qryLenVec = qryLenVecVec[i];
                vector<int> gtVec = gtVecVec[i];
                vector<float> depthVec = depthVecVec[i];

                int64_t indexLeft = -1;
                int64_t indexRight = -1;

                for (size_t j = 0; j < qryLenVec.size(); j++) {
                    uint32_t qryLen = qryLenVec[j];
                    int gt = gtVec[j];
                    float depth = depthVec[j];

                    if (indexLeft == -1 && indexRight == -1) {
                        if (refLen <= 49 && qryLen <= 49) {
                            indexLeft = search_Binary_left(mergeRefStartVec, (refStart-10));
                            indexRight = search_Binary_right(mergeRefStartVec, (refStart+10));
                        } else {
                            indexLeft = search_Binary_left(mergeRefStartVec, (refStart-200));
                            indexRight = search_Binary_right(mergeRefStartVec, (refStart+200));
                        }

                        if (indexLeft < 0 || indexRight >= static_cast<int64_t>(mergeRefStartVec.size())) {
                            throw runtime_error("out of index, please check the data or code.");
                        }
                    }

                    // ------ Match by variant class: breakpoint distance and 25% length difference ------
                    for (int64_t k = indexLeft; k <= indexRight; k++) {
                        vector<uint32_t> mergeQryLenVec = mergeQryLenVecVec[k];
                        uint32_t trueRefLen = mergeRefLenVec[k];
                        uint32_t trueRefStart = mergeRefStartVec[k];

                        vector<uint32_t> trueSeqLenVec;
                        trueSeqLenVec.push_back(trueRefLen);
                        trueSeqLenVec.insert(trueSeqLenVec.end(), mergeQryLenVec.begin(), mergeQryLenVec.end());

                        if (mergeQryLenVec.size() == 0) {
                            throw runtime_error("mergeQryLenVec.size() == 0 -> " + chromosome
                                + " " + to_string(mergeRefStartVec[k]));
                        }

                        for (size_t l = 0; l < trueSeqLenVec.size(); l++) {
                            uint32_t trueQryLen = trueSeqLenVec[l];

                            if (refLen >= 50 && qryLen < 50) {
                                if ((abs(static_cast<int32_t>(refStart) - static_cast<int32_t>(trueRefStart)) <= 200) &&
                                    (abs((static_cast<int32_t>(refStart + refLen) - static_cast<int32_t>(trueRefStart + trueRefLen)) <= 200)) &&
                                    ((abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / (float)trueRefLen) <= 0.25))
                                {
                                    int recall_state;
                                    tie(recall_state, j) = recall_push(
                                        chromosome,
                                        trueRefStart,
                                        software,
                                        qryLenVec,
                                        gtVec,
                                        gt,
                                        depth,
                                        j,
                                        k,
                                        l,
                                        indexLeft,
                                        indexRight
                                    );

                                    if (recall_state == -1)
                                    {
                                        continue;
                                    }

                                    goto stop;
                                }
                            } else if (refLen < 50 && qryLen >= 50) {
                                if ((abs(static_cast<int32_t>(refStart) - static_cast<int32_t>(trueRefStart))<=200) &&
                                    (abs(static_cast<int32_t>(refStart+refLen)-static_cast<int32_t>(trueRefStart+trueRefLen))<=200) &&
                                    ((abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen))/(float)trueQryLen)<=0.25))
                                {
                                    int recall_state;
                                    tie(recall_state, j) = recall_push(
                                        chromosome,
                                        trueRefStart,
                                        software,
                                        qryLenVec,
                                        gtVec,
                                        gt,
                                        depth,
                                        j,
                                        k,
                                        l,
                                        indexLeft,
                                        indexRight
                                    );

                                    if (recall_state == -1)
                                    {
                                        continue;
                                    }

                                    goto stop;
                                }
                            }
                            else if (refLen >= 50 && qryLen >= 50) {
                                if ((abs(static_cast<int32_t>(refStart) - static_cast<int32_t>(trueRefStart))<=200) &&
                                    (abs(static_cast<int32_t>(refStart+refLen)-static_cast<int32_t>(trueRefStart+trueRefLen))<=200) &&
                                    ((abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen))/(float)trueRefLen)<=0.25) &&
                                    ((abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen))/(float)trueQryLen)<=0.25))
                                {
                                    int recall_state;
                                    tie(recall_state, j) = recall_push(
                                        chromosome,
                                        trueRefStart,
                                        software,
                                        qryLenVec,
                                        gtVec,
                                        gt,
                                        depth,
                                        j,
                                        k,
                                        l,
                                        indexLeft,
                                        indexRight
                                    );

                                    if (recall_state == -1)
                                    {
                                        continue;
                                    }

                                    goto stop;
                                }
                            } else if (refLen == 1 && qryLen == 1) {
                                if ((abs(static_cast<int32_t>(refStart) - static_cast<int32_t>(trueRefStart))<=1) &&
                                    (abs(static_cast<int32_t>(refStart+refLen)-static_cast<int32_t>(trueRefStart+trueRefLen))<=1) &&
                                    ((abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen))/(float)trueRefLen)<=0.25) &&
                                    ((abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen))/(float)trueQryLen)<=0.25))
                                {
                                    int recall_state;
                                    tie(recall_state, j) = recall_push(
                                        chromosome,
                                        trueRefStart,
                                        software,
                                        qryLenVec,
                                        gtVec,
                                        gt,
                                        depth,
                                        j,
                                        k,
                                        l,
                                        indexLeft,
                                        indexRight
                                    );

                                    if (recall_state == -1)
                                    {
                                        continue;
                                    }

                                    goto stop;
                                }
                            } else if (refLen >= 3 && qryLen <= 2) {
                                if ((abs(static_cast<int32_t>(refStart) - static_cast<int32_t>(trueRefStart))<=10) &&
                                (abs(static_cast<int32_t>(refStart+refLen)-static_cast<int32_t>(trueRefStart+trueRefLen))<=10) &&
                                ((abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen))/(float)trueRefLen)<=0.25))
                                {
                                    int recall_state;
                                    tie(recall_state, j) = recall_push(
                                        chromosome,
                                        trueRefStart,
                                        software,
                                        qryLenVec,
                                        gtVec,
                                        gt,
                                        depth,
                                        j,
                                        k,
                                        l,
                                        indexLeft,
                                        indexRight
                                    );

                                    if (recall_state == -1)
                                    {
                                        continue;
                                    }

                                    goto stop;
                                }
                            } else if (refLen <= 2 && qryLen >= 3) {
                                if ((abs(static_cast<int32_t>(refStart) - static_cast<int32_t>(trueRefStart))<=10) &&
                                (abs(static_cast<int32_t>(refStart+refLen)-static_cast<int32_t>(trueRefStart+trueRefLen))<=10) &&
                                ((abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen))/(float)trueQryLen)<=0.25))
                                {
                                    int recall_state;
                                    tie(recall_state, j) = recall_push(
                                        chromosome,
                                        trueRefStart,
                                        software,
                                        qryLenVec,
                                        gtVec,
                                        gt,
                                        depth,
                                        j,
                                        k,
                                        l,
                                        indexLeft,
                                        indexRight
                                    );

                                    if (recall_state == -1)
                                    {
                                        continue;
                                    }

                                    goto stop;
                                }
                            } else {
                                if ((abs(static_cast<int32_t>(refStart) - static_cast<int32_t>(trueRefStart))<=1) &&
                                (abs(static_cast<int32_t>(refStart+refLen)-static_cast<int32_t>(trueRefStart+trueRefLen))<=1) &&
                                ((abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen))/(float)trueRefLen)<=0.25) &&
                                ((abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen))/(float)trueQryLen)<=0.25))
                                {
                                    int recall_state;
                                    tie(recall_state, j) = recall_push(
                                        chromosome,
                                        trueRefStart,
                                        software,
                                        qryLenVec,
                                        gtVec,
                                        gt,
                                        depth,
                                        j,
                                        k,
                                        l,
                                        indexLeft,
                                        indexRight
                                    );

                                    if (recall_state == -1)
                                    {
                                        continue;
                                    }

                                    goto stop;
                                }
                            }
                        }
                    }
                    stop:;
                }
            }
        }
    }
}

// Return the length difference between the last non-reference allele of gtVec and REF.
int64_t VCFMerge::sv_length_select(
    const string & vcfInfo,
    const vector<int> & gtVec
) {
    int64_t svLength = 0;

    vector<string> vcfInfoVec = split(vcfInfo, "\t");

    uint32_t refLen;
    vector<uint32_t> qryLenVec;
    tie(refLen, qryLenVec) = get_hap_len(
        vcfInfoVec[2],
        vcfInfoVec[3],
        vcfInfoVec[4],
        gtVec,
        "hap"
    );

    for (size_t i = 0; i < gtVec.size(); i++) {
        if (gtVec[i] == 0) {
            continue;
        }
        else {
            svLength = static_cast<int64_t>(qryLenVec[i] - refLen);
        }
    }

    return svLength;
}

// Choose the tool whose genotype is kept at a site; return an empty string to drop the site.
string VCFMerge::filter_rule(
    vector<vector<int> > & gtAllVec,
    vector<string> & softwareVec,
    vector<float> & depthNorVec,
    const string & vcfInfo
) {
    if (mMode == "all") {
        if (softwareNum_ >= 2 && softwareVec.size() < 2) {
            return "";
        }
    }

    vector<string> bestSoftwareVec{};
    string selectSoftware;
    int64_t svLength{};
    int gtFrequency = 0;

    // ------ Depth closest to the tool mean, most frequent genotype ------
    float selectDepth = 10000.0;
    int selectDepthIdx = 0;
    for (size_t i = 0; i < depthNorVec.size(); i++) {
        if (abs(depthNorVec[i]) < selectDepth) {
            selectDepth = abs(depthNorVec[i]);
            selectDepthIdx = i;
        }
    }

    for (size_t i = 0; i < gtAllVec.size(); i++) {
        svLength = sv_length_select(
            vcfInfo,
            gtAllVec[i]
        );

        if (count(gtAllVec.begin(), gtAllVec.end(), gtAllVec[i]) > gtFrequency) {
            gtFrequency = count(gtAllVec.begin(), gtAllVec.end(), gtAllVec[i]);
        }
    }

    for (size_t i = 0; i < gtAllVec.size(); i++) {
        if (count(gtAllVec.begin(), gtAllVec.end(), gtAllVec[i]) == gtFrequency) {
            bestSoftwareVec.push_back(softwareVec[i]);
        }
    }
    // ------ Selection rule by variant size and mode ------
    if (-49 <= svLength && svLength <= 49) {
        if (mMode == "specific") {
            selectSoftware = (gtFrequency < 2) ? softwareVec[selectDepthIdx] : bestSoftwareVec[0];
        } else {
            if (softwareNum_ == 1) {
                selectSoftware = bestSoftwareVec[0];
            } else if (softwareNum_ == 2 && paragraphInputBool_) {
                selectSoftware = softwareVec[selectDepthIdx];
            } else {
                selectSoftware = (gtFrequency < 2) ? "" : bestSoftwareVec[0];
            }
        }
    } else {
        if (mMode == "specific") {
            selectSoftware = (gtFrequency < 2) ? softwareVec[selectDepthIdx] : bestSoftwareVec[0];
        } else {
            if (softwareNum_ == 1) {
                selectSoftware = bestSoftwareVec[0];
            } else {
                selectSoftware = (gtFrequency < 2) ? softwareVec[selectDepthIdx] : bestSoftwareVec[0];
            }
        }
    }
    return selectSoftware;
}

// Pick one genotype per population site; in mode all also find the population sample most similar to this one and return false if none shares a genotype.
bool VCFMerge::vcf_merge_filter(bool secondBool)
{
    log_stream() << "Filtering.\n";

    uint64_t allPassNum = 0;
    map<uint32_t, uint64_t> sampleIdxNumMap;

    if (secondBool) {
        map<string, map<int, string > >().swap(outChrStartInfoMap_);
    }

    // ------ One genotype per site ------
    for (const auto& [chromosome, startSoftwareInfoMap] : mergeVcfStruct_.recallSoftwareGtDepVecMap) {
        for (const auto& [refStart, SoftwareInfoMap] : startSoftwareInfoMap) {
            string informations;
            vector<string> gtVec;
            tie(informations, gtVec) = mergeVcfStruct_.BaseInfoTupMap[chromosome][refStart];

            string gtTmpTmp = gtVec[mostSimilarSampleIdx_];
            if (secondBool && (gtTmpTmp == "." || gtTmpTmp == ".|." || gtTmpTmp == "./." || gtTmpTmp == "0" || gtTmpTmp == "0|0" || gtTmpTmp == "0/0")) {
                continue;
            }

            vector<vector<int> > gtAllVec;
            vector<string> softwareVec;
            vector<float> depthSumVec;
            vector<float> depthNorVec;
            vector<vector<float> > depthVecVec;

            for (const auto& [software, Info] : SoftwareInfoMap) {
                depthVecVec.push_back(get<0>(Info));

                float meanDepth;
                float sd;
                if (mergeVcfStruct_.depthMap.find(software) != mergeVcfStruct_.depthMap.end()) {
                    meanDepth = get<0>(mergeVcfStruct_.depthMap[software]);
                    sd = get<2>(mergeVcfStruct_.depthMap[software]);
                } else {
                    throw runtime_error("no coverage information -> " + software);
                }
                float depthSumTmp = accumulate(get<0>(Info).begin(), get<0>(Info).end(), 0.0);
                depthSumVec.push_back(depthSumTmp);

                float depthNor = (depthSumTmp - meanDepth) / (float)sd;
                gtAllVec.push_back(get<1>(Info));
                softwareVec.push_back(software);
                depthNorVec.push_back(depthNor);
            }

            string selectSoftware = filter_rule(
                gtAllVec,
                softwareVec,
                depthNorVec,
                informations
            );

            if (selectSoftware.empty()) {
                continue;
            }

            int softwareIdx = distance(softwareVec.begin(), find(softwareVec.begin(), softwareVec.end(), selectSoftware));
            vector<int> selectGtVec = gtAllVec[softwareIdx];
            float selectDepthSumVec = depthSumVec[softwareIdx];
            float selectDepthNor = depthNorVec[softwareIdx];
            vector<float> selectDepthVec = depthVecVec[softwareIdx];

            string selectGt = join(selectGtVec, "/");
            string selectGtTmp = selectGt;
            selectGt += ":" + selectSoftware
                    + ":" + to_string(selectDepthNor).substr(0, 5)
                    + ":" + to_string(get<0>(mergeVcfStruct_.depthMap[selectSoftware])).substr(0, 5)
                    + ":" + to_string(selectDepthSumVec).substr(0, 5);
            outChrStartInfoMap_[chromosome][refStart] = informations + "\t" + selectGt;

            if (mMode == "all") {
                allPassNum++;
                uint32_t sampleIdx = 0;
                replace(selectGtTmp.begin(), selectGtTmp.end(), '/', '|');
                for (const auto& gtTmp : gtVec) {
                    if (gtTmp == "." || gtTmp == ".|." || gtTmp == "./." || gtTmp == "0" || gtTmp == "0|0" || gtTmp == "0/0") {
                        sampleIdx++;
                        continue;
                    }
                    if (selectGtTmp == gtTmp) {
                        auto& emplacedValue = sampleIdxNumMap.emplace(sampleIdx, 0).first->second;
                        emplacedValue++;
                    }
                    sampleIdx++;
                }
            }
        }
    }

    // ------ Population sample sharing the most selected genotypes ------
    if (mMode == "all") {
        if (sampleIdxNumMap.empty()) {
            warning_stream() << "no population sample shares a selected genotype; the high-confidence VCF is not written." << endl;
            return false;
        }
        uint64_t maxValue = sampleIdxNumMap.begin()->second;
        uint32_t maxKey = sampleIdxNumMap.begin()->first;
        for (auto it = sampleIdxNumMap.begin(); it != sampleIdxNumMap.end(); ++it) {
            if (it->second > maxValue) {
                maxValue = it->second;
                maxKey = it->first;
            }
        }
        diagnostic_stream() << fixed << setprecision(2);
        log_stream() << "The sample most similar to the '" << sampleName_ << "' is: '" << mergeVcfStruct_.sampleNameVec[maxKey] << "(" << maxValue / (float)allPassNum << ")'" << endl;
        log_stream() << "If you firmly believe that '" << sampleName_ << "' and '" << mergeVcfStruct_.sampleNameVec[maxKey] << "' belong to closely related varieties, then you can utilize the upcomming high-confidence VCF file." << endl;
        diagnostic_stream() << defaultfloat << setprecision(6);
        mostSimilarSampleIdx_ = maxKey;
    }
    return true;
}

// Return the mean, variance and standard deviation of the depths.
tuple<float, float, float> VCFMerge::cal_var_sd(const vector<float>& data) const
{
    float sum = std::accumulate(std::begin(data), std::end(data), 0.0);
    float mean = sum / data.size();

    float variance = 0.0;
    std::for_each(std::begin(data), std::end(data), [&](const float d) {
        variance += pow(d-mean, 2);
    });
    variance /= data.size();

    float std_deviation = sqrt(variance);

    return make_tuple(mean, variance, std_deviation);
}

// Write the header and the selected records to the output file or stdout.
void VCFMerge::result_save(
)
{
    log_stream() << "Wrote genotyped variants to '" << mOutputFileName << "'\n";

    SAVE SAVEClass(mOutputFileName, threads_);

    stringstream outStream;
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10;
    outStream.str().reserve(CACHE_SIZE);

    outStream << mergeVcfStruct_.headInfo;

    for (const auto& [_, StartInfoMap] : outChrStartInfoMap_) {
        for (const auto& [_, info] : StartInfoMap) {
            outStream << info + "\n";

            if (outStream.tellp() >= CACHE_SIZE) {
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                outStream.str(string());
                outStream.clear();
            }
        }
    }

    if (outStream.tellp() > 0) {
        string outTxt = outStream.str();
        SAVEClass.save(outTxt);
        outStream.str(string());
        outStream.clear();
    }
}

// Split a path into its directory (with the trailing '/') and file name.
void VCFMerge::extract_directory_FileName(const std::string& path, std::string& directory, std::string& fileName) {
    size_t lastSlashPos = path.find_last_of('/');

    if (lastSlashPos != std::string::npos) {
        directory = path.substr(0, lastSlashPos + 1);

        fileName = path.substr(lastSlashPos + 1);
    } else {
        directory = "";
        fileName = path;
    }
}
