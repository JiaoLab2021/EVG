#include "../include/vcf_win.hpp"

#include <regex>
#include <stdexcept>

#include "../include/logger.hpp"

using namespace std;

// Store the window settings, input files, output name, mode and thread count.
VCFWin::VCFWin(
    const uint32_t& windowSize,
    const uint32_t& stepSize,
    const string& chrLenFilename,
    const string& vcfFilename,
    const string& waysFilename,
    const string& outputFileName,
    const string& mode,
    uint32_t threads
) : windowSize_(windowSize), stepSize_(stepSize), chrLenFilename_(chrLenFilename), vcfFilename_(vcfFilename), waysFilename_(waysFilename), outputFileName_(outputFileName), mode_(mode), threads_(threads) {}

// Read the chromosome lengths; throw if the file is missing or a line has fewer than two fields.
void VCFWin::load_chromosome_lengths()
{
    log_stream() << "Count the number and names of chromosomes.\n";

    string line;

    ifstream chrLenFile;
    chrLenFile.open(chrLenFilename_, ios::in);

    if(!chrLenFile.is_open()) {
        throw runtime_error("'" + chrLenFilename_ + "': No such file or directory.");
    } else {
        while(getline(chrLenFile,line)) {
            if (line.empty()) {
                continue;
            } else {
                static const std::regex spaceReg("\\s+");
                std::vector<std::string> lineList(std::sregex_token_iterator(line.begin(), line.end(), spaceReg, -1), std::sregex_token_iterator());

                if (lineList.size() < 2) {
                    throw runtime_error("Format error '" + chrLenFilename_ + "': chrName\tchrLen");
                }

                chrLenMap_[lineList[0]] = stoul(lineList[1]);
            }
        }
    }
    chrLenFile.close();
}

// Build the window starts and ends of every chromosome, plus a last window up to the chromosome end.
void VCFWin::build_windows()
{
    log_stream() << "Building step index.\n";

    for (const auto& [chrName, chrLen] : chrLenMap_) {
        uint32_t stepsNumber = static_cast<uint32_t>(std::ceil(chrLen / stepSize_));

        for (uint32_t i = 0; i < stepsNumber; i++) {
            uint32_t stepStart = i * stepSize_ + 1;
            uint32_t stepEnd = stepStart + windowSize_ - 1;
            winStartMap_[chrName].push_back(stepStart);
            winEndMap_[chrName].push_back(stepEnd);
        }
        if (winEndMap_[chrName].back()<chrLen) {
            winStartMap_[chrName].push_back(winEndMap_[chrName].back() + 1);
            winEndMap_[chrName].push_back(chrLen);
        }
    }
}

// Read the sample groups; without -g every sample belongs to group 0.
void VCFWin::load_groups()
{
    log_stream() << "Building group index.\n";

    if (waysFilename_.empty()) {
        waysGroupMap_["0"] = "0";
        return;
    }

    ifstream waysFile;
    waysFile.open(waysFilename_, ios::in);

    string line;
    if (!waysFile.is_open()) {
        throw runtime_error("'" + waysFilename_ + "': No such file or directory.");
    } else {
        while (getline(waysFile, line)) {
            if (line.empty()) {
                continue;
            } else {
                line = strip(line, '\n');

                static const std::regex spaceReg("\\s+");
                std::vector<std::string> lineList(std::sregex_token_iterator(line.begin(), line.end(), spaceReg, -1), std::sregex_token_iterator());

                if (lineList.size() < 2) {
                    lineList = split(line, " ");
                }

                if (lineList.size() < 2) {
                    throw runtime_error("Format error '" + waysFilename_ + "': way\tgroup");
                }

                waysGroupMap_[lineList[0]] = lineList[1];
            }
        }
    }

    waysFile.close();
}

// Add the length and group of every non-reference sample genotype to the windows overlapping the variant start.
void VCFWin::count_window_variants()
{
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFilename_, threads_);

    vector<string> waysVector;

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#CHROM") != string::npos) {
            waysVector = INFOSTRUCTTMP.lineVec;
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        if (waysVector.size() < 1) {
            throw runtime_error("missing the comment line: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  CW01    CW02......");
        }

        vector<uint32_t> qryLenList;
        for (size_t i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++) {
            qryLenList.push_back(INFOSTRUCTTMP.ALTVec[i].size());
        }

        for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
            vector<int> gtVec = get_gt(
                INFOSTRUCTTMP.lineVec,
                i
            );

            string gt = join(gtVec, "/");

            if (gt == "." || gt == "0" || gt == "0/0" || gt =="0|0") {
                continue;
            }

            vector<string> stepsVectorKey;

            for (auto rit = winStartMap_.begin(); rit != winStartMap_.end(); ++rit) {
                stepsVectorKey.push_back(rit->first);
            }

            if (find(stepsVectorKey.begin(), stepsVectorKey.end(), INFOSTRUCTTMP.CHROM) == stepsVectorKey.end()) {
                throw runtime_error(INFOSTRUCTTMP.CHROM + "\t" + (stepsVectorKey.empty() ? string() : stepsVectorKey[0]));
            } else {
                uint32_t indexRight = search_Binary_right(winEndMap_[INFOSTRUCTTMP.CHROM], INFOSTRUCTTMP.POS);
                uint32_t indexLeft = search_Binary_left(winStartMap_[INFOSTRUCTTMP.CHROM], INFOSTRUCTTMP.POS);

                string way = waysVector[i];
                string group;
                if (waysGroupMap_.find(way) == waysGroupMap_.end()) {
                    group = "0";
                } else {
                    group = waysGroupMap_[way];
                }

                uint32_t winStart1 = winStartMap_[INFOSTRUCTTMP.CHROM][indexRight];
                uint32_t winStart2 = winStartMap_[INFOSTRUCTTMP.CHROM][indexLeft];

                if (indexLeft - indexRight == 1) {
                    chrWinStartLenMap_[INFOSTRUCTTMP.CHROM][winStart1].push_back(min(winEndMap_[INFOSTRUCTTMP.CHROM][indexRight], INFOSTRUCTTMP.END)-INFOSTRUCTTMP.POS+1);
                    chrWinStartGroupMap_[INFOSTRUCTTMP.CHROM][winStart1].push_back(group);
                }

                chrWinStartLenMap_[INFOSTRUCTTMP.CHROM][winStart2].push_back(min(winStart2+windowSize_-1, INFOSTRUCTTMP.END)-INFOSTRUCTTMP.POS+1);
                chrWinStartGroupMap_[INFOSTRUCTTMP.CHROM][winStart2].push_back(group);
            }
        }
    }
}

// Write the per-window sum, average and per-group fractions.
void VCFWin::save_result()
{
    int waysNumber = waysGroupMap_.size();

    map<string, int> groupNumMap;
    for (auto rit = waysGroupMap_.begin(); rit != waysGroupMap_.end(); ++rit) {
        string group = rit->second;
        groupNumMap[group]++;
    }

    SAVE outFile(outputFileName_, threads_);

    string ouTxt = "#waysNumber: " + to_string(waysNumber) + "\n";
    ouTxt += "#chr\tLocation\t" + mode_ + "\taverage";
    for (auto iter1 : groupNumMap) {
        ouTxt += "\t" + iter1.first;
    }
    ouTxt += "\n";

    long double sumValue;
    long double averageValue;

    for (const auto& [chromosome, winStartVec] : winStartMap_) {
        int idxTmp = 0;

        if (chrWinStartLenMap_.find(chromosome) == chrWinStartLenMap_.end()) {
            continue;
        }

        for (const auto& winStart : winStartVec) {
            uint32_t winEnd = winEndMap_[chromosome][idxTmp];
            map<uint32_t, vector<uint32_t> >::iterator findIter1 = chrWinStartLenMap_[chromosome].find(winStart);
            if (findIter1 != chrWinStartLenMap_[chromosome].end()) {
                if (mode_ == "length") {
                    sumValue = accumulate(findIter1->second.begin(), findIter1->second.end(), 0);
                    averageValue = sumValue/(winEnd-winStart+1);
                } else {
                    sumValue = findIter1->second.size();
                    averageValue = sumValue/waysNumber;
                }

                ouTxt += chromosome + "\t" + to_string(winStart) + "-" + to_string(winEnd) + "\t" + to_string(sumValue) + "\t" + to_string(averageValue);

                for (auto iter3 : groupNumMap) {
                    string group = iter3.first;
                    int groupNum = iter3.second;
                    int number = count(chrWinStartGroupMap_[chromosome][winStart].begin(), chrWinStartGroupMap_[chromosome][winStart].end(), group);
                    ouTxt += "\t" + to_string(number/(float)groupNum);
                }
            } else {
                ouTxt += chromosome + "\t" + to_string(winStart) + "-" + to_string(winEnd) + "\t0\t0";

                for (auto iter3 : groupNumMap) {
                    ouTxt += "\t0";
                }
            }
            ouTxt += "\n";
            idxTmp++;
        }
    }
    outFile.save(ouTxt);
}

// Return the genotype of one sample (0 = last column) as integers; missing or unseparated genotypes become {0, 0}.
vector<int> VCFWin::get_gt(
    const vector<string> & lineVec,
    int sampleIdx
)
{
    vector <int> gtVec;

    int formatIndex = 8;

    vector<string> formatVec = split(lineVec[formatIndex], ":");

    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));

    if (gtIndex == formatVec.size()) {
        warning_stream() << "no [GT] information in FORMAT -> " << lineVec[0] << ":" << lineVec[1] << endl;
        gtVec = {0, 0};
    } else {
        string gt;

        if (sampleIdx == 0) {
            gt = VCFOPEN::sample_gt(lineVec, lineVec.size() - 1, gtIndex);
        } else {
            gt = VCFOPEN::sample_gt(lineVec, sampleIdx, gtIndex);
        }

        string splitStr;
        if (gt.find("/") != string::npos) {
            splitStr = "/";
        } else if (gt.find("|") != string::npos) {
            splitStr = "|";
        } else {
            gtVec = {0, 0};
            return gtVec;
        }

        for (auto it : split(gt, splitStr)) {
            if (it == ".") {
                gtVec = {0, 0};
                return gtVec;
            }
            gtVec.push_back(stoi(it));
        }
    }

    return gtVec;
}
