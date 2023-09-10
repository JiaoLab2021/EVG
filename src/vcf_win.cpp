// g++ vcf_win_count_thread.cpp -std=c++11 -o vcf_win_count_thread -lpthread -O2
#include "../include/vcf_win.hpp"

using namespace std;

int main_win(int argc, char** argv)
{
    string vcfFilename;
    string chrLenFilename;
    string waysFilename;

    string outputFileName;

    uint32_t windowSize = 200000;
    uint32_t stepSize = 100000;
    string mode = "number";

    // Input parameter
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},
            {"length", required_argument, 0, 'l'},         
            {"group", required_argument, 0, 'g'},
            {"out", required_argument, 0, 'o'},

            {"windowSize", required_argument, 0, 'w'},
            {"stepSize", required_argument, 0, 's'},
            {"mode", required_argument, 0, 'm'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:l:g:o:w:s:m:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFilename = optarg;
            break;
        case 'l':
            chrLenFilename = optarg;
            break;
        case 'g':
            waysFilename = optarg;
            break;
        case 'o':
            outputFileName = optarg;
            break;
        case 'w':
            windowSize = stoi(optarg);
            break;
        case 's':
            stepSize = stoi(optarg);
            break;
        case 'm':
            mode = optarg;
            break;
        case 'h':
        case '?':
            help_win(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 2) {
        help_win(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // init
    VCFWin VCFWinClass(
        windowSize, 
        stepSize, 
        chrLenFilename, 
        vcfFilename, 
        waysFilename, 
        outputFileName,  
        mode
    );

    // Count chromosome numbers and names
    VCFWinClass.count_chrName();

    // Step dictionary
    VCFWinClass.step_count();

    // Classification of ways
    VCFWinClass.ways_group();
    
    // Calculate the length of the variation in the window
    VCFWinClass.window_len_count();

    // save
    VCFWinClass.save_result();

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";
    return 0;
}

// Help document
void help_win(char** argv)

{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v -l [options]" << endl
       << "calculate the variants density in VCF file using window-based statistics." << endl
       << endl
       << "input/output:" << endl
       << "    -v, --vcf         FILE        vcf file merged by vcftools" << endl
       << "    -l, --length      FILE        chromosome length file (chr\tlength)" << endl
       << "    -o, --out         FILE        output genotyping to FILE [stdout]" << endl
       << endl
       << "optional parameter:" << endl
       << "    -w, --windowSize  INT         window [200000]" << endl
       << "    -s, --stepSize    INT         step [100000]" << endl
       << "    -g, --group       FILE        group information for every way (way\tgroup)" << endl
       << "    -m, --mode        STRING      which indicator to count (length/number) [number]" << endl
       << endl
       << "    -h, --help                    print this help document" << endl;
}


/**
 * init
 *
 * @param windowSize
 * @param stepSize
 * @param chrLenFilename
 * @param vcfFilename
 * @param waysFilename
 * @param outputFileName
 * @param mode
 * 
**/
VCFWin::VCFWin(
    const uint32_t& windowSize, 
    const uint32_t& stepSize, 
    const string& chrLenFilename, 
    const string& vcfFilename, 
    const string& waysFilename, 
    const string& outputFileName, 
    const string& mode
) : windowSize_(windowSize), stepSize_(stepSize), chrLenFilename_(chrLenFilename), vcfFilename_(vcfFilename), waysFilename_(waysFilename), outputFileName_(outputFileName), mode_(mode) {}


// Count chromosome numbers and names
void VCFWin::count_chrName()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Count the number and names of chromosomes.\n";

    // See how many chromosomes there are
    string line;

    // input file stream
    ifstream chrLenFile;
    chrLenFile.open(chrLenFilename_, ios::in);

    if(!chrLenFile.is_open()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << chrLenFilename_ << "': No such file or directory." << endl;
        exit(1);
    } else {
        while(getline(chrLenFile,line)) {
            if (line.empty()) {
                continue;
            } else {
                std::regex regex("\\s+"); // Use the regular expression "\\s+" to match whitespace characters as delimiters
                std::vector<std::string> lineList(std::sregex_token_iterator(line.begin(), line.end(), regex, -1), std::sregex_token_iterator());

                if (lineList.size() < 2) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " 
                        << "Format error \'" << chrLenFilename_ << "\'" << ": chrName\tchrLen\n";
                    exit(1);
                }
                
                chrLenMap_[lineList[0]] = stoul(lineList[1]);    
            }
        }
    }
    // Free memory
    chrLenFile.close();
}

// Step dictionary
void VCFWin::step_count()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building step index.\n";

    for (const auto& [chrName, chrLen] : chrLenMap_) {
        // Step number
        uint32_t stepsNumber = static_cast<uint32_t>(std::ceil(chrLen / stepSize_));

        // Make a step dictionary
        for (uint32_t i = 0; i < stepsNumber; i++) {
            uint32_t stepStart = i * stepSize_ + 1;
            uint32_t stepEnd = stepStart + windowSize_ - 1;
            winStartMap_[chrName].push_back(stepStart);
            winEndMap_[chrName].push_back(stepEnd);
        }
        // And finally to the end of the chromosome
        if (winEndMap_[chrName].back()<chrLen) {
            winStartMap_[chrName].push_back(winEndMap_[chrName].back() + 1);
            winEndMap_[chrName].push_back(chrLen);
        }
    }
}

// Classification of ways
void VCFWin::ways_group()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building group index.\n";

    if (waysFilename_.empty()) {  // If the group file is empty, a fake group hash table is built
    
        waysGroupMap_["0"] = "0";
        return;
    }

    ifstream waysFile;
    waysFile.open(waysFilename_, ios::in);

    string line;
    if (!waysFile.is_open()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << waysFilename_ << "': No such file or directory." << endl;
        exit(1);
    } else {
        while (getline(waysFile, line)) {
            if (line.empty()) {
                continue;
            } else {
                line = strip(line, '\n');

                std::regex regex("\\s+"); // Use the regular expression "\\s+" to match whitespace characters as delimiters
                std::vector<std::string> lineList(std::sregex_token_iterator(line.begin(), line.end(), regex, -1), std::sregex_token_iterator());

                if (lineList.size() < 2) {
                    lineList = split(line, " ");
                }

                if (lineList.size() < 2) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Format error \'" << waysFilename_ << "\': way\tgroup\n";
                    exit(1);
                }

                waysGroupMap_[lineList[0]] = lineList[1];
            }
        }
    }

    // Free memory
    waysFile.close();
}


// Calculate the length of the variation in the window
void VCFWin::window_len_count()
{
     // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFilename_);

    vector<string> waysVector;

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // skip empty line
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
        
        // If the vcf file has no comment line, the code simply exits.
        if (waysVector.size() < 1) {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: missing the comment line: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  CW01    CW02......\n";
            exit(1);
        }

        // Split the qry sequence by ',' and loop it into qryLenList
        vector<uint32_t> qryLenList;
        for (size_t i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++) {
            qryLenList.push_back(INFOSTRUCTTMP.ALTVec[i].size());
        }

        // The vcf is genotype information starting from column 10, so the loop starts from i=9, which is column 10
        // Different vcf file formats may be different and need to be checked
        for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
            vector<int> gtVec = get_gt(
                INFOSTRUCTTMP.lineVec, 
                i
            );

            string gt = join(gtVec, "/");

            // Skip if gt is empty
            if (gt == "." || gt == "0" || gt == "0/0" || gt =="0|0") {
                continue;
            }

            vector<string> stepsVectorKey;

            // In the index of the step dictionary, check whether the chromosome is in the dictionary first, and report an error if it is not
            for (auto rit = winStartMap_.begin(); rit != winStartMap_.end(); ++rit) {
                stepsVectorKey.push_back(rit->first);                       
            }

            if (find(stepsVectorKey.begin(), stepsVectorKey.end(), INFOSTRUCTTMP.CHROM) == stepsVectorKey.end()) {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "IndexError: " << INFOSTRUCTTMP.CHROM << "\t" << stepsVectorKey[0] << endl;
                exit(1);
            } else {
                uint32_t indexRight = search_Binary_right(winEndMap_[INFOSTRUCTTMP.CHROM], INFOSTRUCTTMP.POS);  // Previous window coordinates
                uint32_t indexLeft = search_Binary_left(winStartMap_[INFOSTRUCTTMP.CHROM], INFOSTRUCTTMP.POS);  // Next window coordinates
                
                string way = waysVector[i];  // The way of the strain
                string group;  // The group of the strain
                if (waysGroupMap_.find(way) == waysGroupMap_.end()) {  // group
                    group = "0";
                } else {
                    group = waysGroupMap_[way];
                }

                uint32_t winStart1 = winStartMap_[INFOSTRUCTTMP.CHROM][indexRight];  // Previous window
                uint32_t winStart2 = winStartMap_[INFOSTRUCTTMP.CHROM][indexLeft];  // This window

                // This window
                if (indexLeft - indexRight == 1) {  // Previous window info added
                    chrWinStartLenMap_[INFOSTRUCTTMP.CHROM][winStart1].push_back(min(winEndMap_[INFOSTRUCTTMP.CHROM][indexRight], INFOSTRUCTTMP.END)-INFOSTRUCTTMP.POS+1);  // length
                    chrWinStartGroupMap_[INFOSTRUCTTMP.CHROM][winStart1].push_back(group);  // group
                }

                // ±¾´°¿Ú
                chrWinStartLenMap_[INFOSTRUCTTMP.CHROM][winStart2].push_back(min(winStart2+windowSize_-1, INFOSTRUCTTMP.END)-INFOSTRUCTTMP.POS+1);  // length
                chrWinStartGroupMap_[INFOSTRUCTTMP.CHROM][winStart2].push_back(group);  // group
            }
        }
    }
}


/**
 * save result
 * 
 * 
 * @return void
**/
void VCFWin::save_result()
{
    int waysNumber = waysGroupMap_.size();

    // Record the number of each group <group, number>
    map<string, int> groupNumMap;
    for (auto rit = waysGroupMap_.begin(); rit != waysGroupMap_.end(); ++rit) {
        string group = rit->second;
        groupNumMap[group]++;
    }

    // save result
    SAVE outFile(outputFileName_);

    // Print head
    string ouTxt = "#waysNumber: " + to_string(waysNumber) + "\n";
    ouTxt += "#chr\tLocation\t" + mode_ + "\taverage";
    for (auto iter1 : groupNumMap) {
        ouTxt += "\t" + iter1.first;
    }
    ouTxt += "\n";

    // Calculation result
    long double sumValue;
    long double averageValue;

    for (const auto& [chromosome, winStartVec] : winStartMap_) {  // map<string, vector<start> >
        int idxTmp = 0;

        // Threads do not print if their chromosomes do not correspond
        if (chrWinStartLenMap_.find(chromosome) == chrWinStartLenMap_.end()) {
            continue;
        }

        for (const auto& winStart : winStartVec) {
            uint32_t winEnd = winEndMap_[chromosome][idxTmp];
            map<uint32_t, vector<uint32_t> >::iterator findIter1 = chrWinStartLenMap_[chromosome].find(winStart);
            if (findIter1 != chrWinStartLenMap_[chromosome].end()) {
                // The length of the variation in the statistical window
                if (mode_ == "length") {
                    // The sum of the window variation lengths
                    sumValue = accumulate(findIter1->second.begin(), findIter1->second.end(), 0);
                    // The proportion of variation length to window
                    averageValue = sumValue/(winEnd-winStart+1);
                } else {  // Count the number of variations in the window
                    // The sum of the number of variations in the window
                    sumValue = findIter1->second.size();
                    // The average number of variations
                    averageValue = sumValue/waysNumber;
                }
                
                ouTxt += chromosome + "\t" + to_string(winStart) + "-" + to_string(winEnd) + "\t" + to_string(sumValue) + "\t" + to_string(averageValue);
                    
                // Indicates the number of groups
                for (auto iter3 : groupNumMap) {
                    string group = iter3.first;
                    int groupNum = iter3.second;
                    int number = count(chrWinStartGroupMap_[chromosome][winStart].begin(), chrWinStartGroupMap_[chromosome][winStart].end(), group);
                    ouTxt += "\t" + to_string(number/(float)groupNum);
                }
            } else {
                ouTxt += chromosome + "\t" + to_string(winStart) + "-" + to_string(winEnd) + "\t0\t0";

                // Indicates the number of groups
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

/**
 * Indicates the number of groups.
 *
 * @param lineVec          lineVec
 * @param sampleIdx        Index of the sample genotype, with the default value 0 representing the last column
 * 
 * 
 * @return gtVec           vector <int>
**/
vector<int> VCFWin::get_gt(
    const vector<string> & lineVec, 
    int sampleIdx
)
{
    vector <int> gtVec;  // Locus typing vector

    int formatIndex = 8; // FORMAT column

    // FORMAT character split
    vector<string> formatVec = split(lineVec[formatIndex], ":");

    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // Gets the index location of GT

    if (gtIndex == formatVec.size()) {  // Determine whether the index exists. If it does not exist, the genotypes are all 0.
    
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << lineVec[0] << ":" << lineVec[1] << endl;
        gtVec = {0, 0};
    } else {  // If it exists, save it
        string gt;  // Store the genotype field

        if (sampleIdx == 0) {  // The number of columns not specified is the last column
            gt = split(lineVec.back(), ":")[gtIndex];  // gt field
        } else {
            gt = split(lineVec[sampleIdx], ":")[gtIndex];  // gt field
        }
        
        string splitStr;  // Delimiter in gt
        if (gt.find("/") != string::npos) {  // Determine the '/' separator
            splitStr = "/";
        } else if (gt.find("|") != string::npos) {  // Determine that '|' is the separator
            splitStr = "|";
        } else {  // Return a null value if you do not know
        
            gtVec = {0, 0};
            return gtVec;
        }
        
        for (auto it : split(gt, splitStr)) {  // Once you find gt, split it by splitStr and loop
            if (it == ".") {  // If it is '.', skip this site
            
                gtVec = {0, 0};
                return gtVec;
            }
            gtVec.push_back(stoi(it));  // Add to vector
        }
    }

    return gtVec;
}