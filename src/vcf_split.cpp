// g++ vcf_split.cpp -o vcf_split -lz -O2
#include <getopt.h>
#include <cstdlib>

#include "../include/vcf_split.hpp"


using namespace std;

int main_split(int argc, char** argv)
{
    // Input file
    string vcfFileName;

    // Input to extract the vcf file for snp and indel locations
    string baseVcfFileName;
    
    // Split mode
    string mode = "type";

    // The distance of snp and indel from sv
    int length = 100;

    // Output file prefix
    string prefix = "split";

    // Whether to FILTER based on the filter column
    bool filterBool = false;

    // Input parameter
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"mode", required_argument, 0, 'm'},
            {"length", required_argument, 0, 'l'},
            {"basevcf", required_argument, 0, 'V'},
            {"prefix", required_argument, 0, 'p'},
            {"filter", no_argument, 0, 'f'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:m:l:V:p:fh", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'm':
            mode = optarg;
            break;
        case 'l':
            length = stoi(optarg);
            break;
        case 'V':
            baseVcfFileName = optarg;
            break;
        case 'p':
            prefix = optarg;
            break;;
        case 'f':
            filterBool = true;
            break;
        case 'h':
        case '?':
            help_split(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // Check whether the parameters are correct
    if (vcfFileName.empty() || (mode != "type" && mode != "number"))
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_split(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // init
    VCFSplit VCFSplitClass(vcfFileName, baseVcfFileName, length, filterBool, prefix);
    
    if (mode == "type")
    {
        VCFSplitClass.vcf_split_type();
    }
    else if (mode == "number")
    {
        VCFSplitClass.vcf_split_number();
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_split(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// Help document
void help_split(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "split VCF files by type or number of nearby SNP+Indels." << endl
       << endl
       << "required arguments (Note: vcf files must be sorted):" << endl
       << "    -v, --vcf           FILE       vcf file to be converted" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -m, --mode          STRING     split standard: (type -> SNP+Indel+Ins+Del / number -> number of nearby SNP+Indel) [type]" << endl
       << "    -l, --length        INT        distance of snp and indel from sv [100]" << endl      
       << "    -V, --basevcf       INT        vcf file used to extract snp and indel locations [same as -v]" << endl      
       << "    -p, --prefix        STRING     output prefix [split]" << endl
       << "    -f, --filter        BOOL       filter using the FILTER column" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}



/**
 * @brief Convert vcf to the format required by graph genome tools
 * 
 * @param vcfFileName      input VCF file name
 * @param baseVcfFileName  VCF file for extracting snp and indel locations
 * @param length           Count the number of SNPs and indels within the 'length' near the breakpoint
 * @param filterBool       Whether to filter according to the 'FILTER' column
 * @param prefix           output file prefix
 * 
**/
VCFSplit::VCFSplit(
    string vcfFileName,
    string baseVcfFileName,
    int length, 
    bool filterBool,
    string prefix
) : vcfFileName_(vcfFileName), baseVcfFileName_(baseVcfFileName), length_(length), filterBool_(filterBool), prefix_(prefix) {}


/**
 * @brief Get the length information of variant
 * 
 * @param INFOSTRUCTTMP   variant information
 * @param VCFOPENCLASS    VCFOPEN
 * @param refLen          store the length of REF
 * @param qryLen          store the length of ALT
 * 
**/
bool VCFSplit::get_len(
    const VCFINFOSTRUCT& INFOSTRUCTTMP, 
    VCFOPEN& VCFOPENCLASS, 
    uint32_t& refLen, 
    uint32_t& qryLen
)
{
    // The last column of information
    vector<string> lastColVec = split(INFOSTRUCTTMP.lineVec[INFOSTRUCTTMP.lineVec.size()-1], ":");

    // Find the location of gt in the FORMAT field
    vector<string> formatVec = split(strip(INFOSTRUCTTMP.lineVec[8], '\n'), ":");
    vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");

    string gt;
    int gtIndex = 0;
    if (gtItera != formatVec.end()) {  // The FORMAT contains GT
        // GT index
        gtIndex = distance(formatVec.begin(), gtItera);
        gt = strip(lastColVec[gtIndex], '\n');
    } else {  // If not, exit the code
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "Error: GT not in FORMAT column -> " 
            << INFOSTRUCTTMP.line << endl;
        exit(1);
    }

    if (gt == "0/0" || gt == "." || gt == "./." || gt == ".|." || gt == "0|0") {  // If it is (0/0,.) Format is directly skipped, do not count.
        return false;
    }

    // According to the FILTER field
    if (INFOSTRUCTTMP.FILTER != "PASS" && filterBool_) {  // If it is not (PASS), it is skipped. If it is not filtered, it is skipped
        return false;
    }

    // Software classification results
    vector<string> evaluate_gt;
    if (filterBool_) {
        evaluate_gt = VCFOPENCLASS.gt_split(gt);
    }
    else {  // If not filtered, then there is no judgment here either, because gramtools' GT format is 1 with no separator
        evaluate_gt = {"1", "1"};
    }
    
    if (INFOSTRUCTTMP.ALT.find(">") != string::npos) {  // The result of graphtyper software
        // Use regular expressions to extract length information
        std::regex reg("SVSIZE=(\\d+)");
        std::smatch match;
        // Check whether there is any length information in the qry sequence. If there is no length information, skip this site
        if (!std::regex_search(INFOSTRUCTTMP.ALT, match, reg)) {  // <INS:SVSIZE=97:BREAKPOINT1>
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: No length information in the fourth column, skip this site -> " << INFOSTRUCTTMP.line << endl; 
            return false;
        }

        if (INFOSTRUCTTMP.ALT.find("<INS") != string::npos) {  // Character segment containing inserted fields (graphtyper) <INS:SVSIZE=97:BREAKPOINT1>
            refLen = INFOSTRUCTTMP.LEN;
            qryLen = std::stoul(match[1].str());
        } else if (INFOSTRUCTTMP.ALT.find("<DEL") != string::npos) {  // Character segment contains deletion fields (graphtyper) <DEL:SVSIZE=233:COVERAGE>
            refLen = std::stoul(match[1].str());
            qryLen = INFOSTRUCTTMP.LEN;
        } else if (INFOSTRUCTTMP.ALT.find("<DUP") != string::npos) {  // Character segment containing duplicate fields (graphtyper) <DUP:SVSIZE=2806:BREAKPOINT1>
            refLen = std::stoul(match[1].str());
            qryLen = refLen * 2;
        } else {  // Unknown fields
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: unknown string -> ref_seq:" << INFOSTRUCTTMP.REF << " qey_seq:" << INFOSTRUCTTMP.ALT << endl;
            refLen = INFOSTRUCTTMP.LEN;

            qryLen = std::stoul(match[1].str());
        }
    } else {  // Normal variation
        refLen = INFOSTRUCTTMP.LEN;

        string qrySeq;
        
        uint32_t maxGtNum = stoi(*max_element(evaluate_gt.begin(), evaluate_gt.end()));
        if (maxGtNum > INFOSTRUCTTMP.ALTVec.size()) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: genotyping and ALT sequence numbers do not match -> " << INFOSTRUCTTMP.line << endl;
            exit(1);
        } else {
            qrySeq = INFOSTRUCTTMP.ALTVec[maxGtNum-1]; // GT number minus 1 is the index of the sequence
        }

        qryLen = qrySeq.size();
    }

    // The result of bayestyper was judged to be duplication, ref_len was 1-2, and qry_seq was very long
    if ((INFOSTRUCTTMP.ID.find("Duplication") != string::npos) && (refLen <= 2)) {  // bayestyper will change duplication into insertion, so recalculate the length
        refLen = qryLen;
        qryLen *= 2;
    }

    return true;
}


void VCFSplit::vcf_split_type()
{
    // Input file stream
    // Store vcf information
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    stringstream simpleSVsStream, invSVsStream, dupSVsStream, otherSVsStream; // Use stringstream instead of string concatenation
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
    simpleSVsStream.str().reserve(CACHE_SIZE);
    invSVsStream.str().reserve(CACHE_SIZE);
    dupSVsStream.str().reserve(CACHE_SIZE);
    otherSVsStream.str().reserve(CACHE_SIZE);

    // Output file stream
    SAVE SimpleSave(prefix_ + ".snp.indel.ins.del.vcf.gz");
    SAVE InvSave(prefix_ + ".inv.vcf.gz");
    SAVE DupSave(prefix_ + ".dup.vcf.gz");
    SAVE OtherSave(prefix_ + ".other.vcf.gz");

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#") != string::npos) {  // Check for comment lines
            string headLine = INFOSTRUCTTMP.line + "\n";
            SimpleSave.save(headLine);
            InvSave.save(headLine);
            DupSave.save(headLine);
            OtherSave.save(headLine);
        } else {
            // Get the length information of variant. If encountering an unrecognized string, continue to the next loop.
            uint32_t refLen, qryLen;
            if (!get_len(INFOSTRUCTTMP, VCFOPENCLASS, refLen, qryLen)) {
                continue;
            }

            // Determine the variation type and save it to the string
            if ((refLen < 50 && qryLen < 50) || (refLen >= 50 && qryLen < 50) || (refLen < 50 && qryLen >= 50)) {  // snp+indel+del+ins
                simpleSVsStream << INFOSTRUCTTMP.line + "\n";
            } else if (refLen >= 50 && qryLen >= 50) {
                if (refLen*2 <= qryLen+10 && refLen*2 >= qryLen-10) {  // dup
                    dupSVsStream << INFOSTRUCTTMP.line + "\n";
                } else if (refLen <= qryLen+10 && refLen >= qryLen-10) {  // inv
                    invSVsStream << INFOSTRUCTTMP.line + "\n";
                } else {
                    otherSVsStream << INFOSTRUCTTMP.line + "\n";
                }
            } else {
                otherSVsStream << INFOSTRUCTTMP.line + "\n";
            }

            // save result
            if (simpleSVsStream.tellp() >= CACHE_SIZE || 
                invSVsStream.tellp() >= CACHE_SIZE || 
                dupSVsStream.tellp() >= CACHE_SIZE || 
                otherSVsStream.tellp() >= CACHE_SIZE
            ) {  // Determine the variation type and save it to the string
                string simpleSVs = simpleSVsStream.str();
                string invSVs = invSVsStream.str();
                string dupSVs = dupSVsStream.str();
                string otherSVs = otherSVsStream.str();

                SimpleSave.save(simpleSVs);
                InvSave.save(invSVs);
                DupSave.save(dupSVs);
                OtherSave.save(otherSVs);

                // Clear stringstream
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

    // Last write
    if (simpleSVsStream.tellp() > 0 || 
        invSVsStream.tellp() > 0 || 
        dupSVsStream.tellp() > 0 || 
        otherSVsStream.tellp() > 0
    ) {  // Determine the variation type and save it to the string
        string simpleSVs = simpleSVsStream.str();
        string invSVs = invSVsStream.str();
        string dupSVs = dupSVsStream.str();
        string otherSVs = otherSVsStream.str();

        SimpleSave.save(simpleSVs);
        InvSave.save(invSVs);
        DupSave.save(dupSVs);
        OtherSave.save(otherSVs);

        // Clear stringstream
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


// Classified by the number of nearby SNPs and indels
void VCFSplit::vcf_split_number()
{
    if (baseVcfFileName_.empty()) {
        baseVcfFileName_ = vcfFileName_;
    }

    stringstream SVs0Stream, SVs1Stream, SVs2Stream, SVs4Stream, SVs6Stream, SVs8Stream, SVsMoreStream; // Use stringstream instead of string concatenation
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
    SVs0Stream.str().reserve(CACHE_SIZE);
    SVs1Stream.str().reserve(CACHE_SIZE);
    SVs2Stream.str().reserve(CACHE_SIZE);
    SVs4Stream.str().reserve(CACHE_SIZE);
    SVs6Stream.str().reserve(CACHE_SIZE);
    SVs8Stream.str().reserve(CACHE_SIZE);
    SVsMoreStream.str().reserve(CACHE_SIZE);

    // Output file stream
    SAVE SVs0Save(prefix_ + ".0.vcf.gz");
    SAVE SVs1Save(prefix_ + ".1.vcf.gz");
    SAVE SVs2Save(prefix_ + ".2.vcf.gz");
    SAVE SVs4Save(prefix_ + ".4.vcf.gz");
    SAVE SVs6Save(prefix_ + ".6.vcf.gz");
    SAVE SVs8Save(prefix_ + ".8.vcf.gz");
    SAVE SVsMoreSave(prefix_ + ".more.vcf.gz");

    // snp+indel index
    VCFINFOSTRUCT INFOSTRUCTBase;
    VCFOPEN VCFOPENCLASSBase(baseVcfFileName_);

    // Record the location information of snp and indel
    unordered_map<string, snpIndelInfo> snpIndelInfoMap;  // unordered_map<chromosome, snpIndelInfo>

    // If not traversed, continue
    while (VCFOPENCLASSBase.read(INFOSTRUCTBase)) {
        // empty line, skip
        if (INFOSTRUCTBase.line.empty()) {
            continue;
        }

        if (INFOSTRUCTBase.line.find("#") == string::npos) {
            // Get the length information of variant. If encountering an unrecognized string, continue to the next loop.
            uint32_t refLen, qryLen;
            if (!get_len(INFOSTRUCTBase, VCFOPENCLASSBase, refLen, qryLen)) {
                continue;
            }

            uint32_t refEnd;
            refEnd = INFOSTRUCTBase.POS + refLen - 1;

            // Determine the variation type and save it to the hash table
            if ((refLen < 50 && qryLen < 50)) {  // Save the start and end locations of snp+indel
                snpIndelInfoMap[INFOSTRUCTBase.CHROM].refStartVec.push_back(INFOSTRUCTBase.POS);
                snpIndelInfoMap[INFOSTRUCTBase.CHROM].refEndVec.push_back(refEnd);
            }
        }
    }


    // Input file stream
    // Store vcf information
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // Building a hash table (structure variation)
    map<string, map<uint32_t, svInfo> > svInfoMap; // map<chromosome, map<refStart, svInfo>>

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#") != string::npos) {  // Check for comment lines
            string headLine = INFOSTRUCTTMP.line + "\n";
            SVs0Save.save(headLine);
            SVs1Save.save(headLine);
            SVs2Save.save(headLine);
            SVs4Save.save(headLine);
            SVs6Save.save(headLine);
            SVs8Save.save(headLine);
            SVsMoreSave.save(headLine);
        } else {
            // Get the length information of variant. If encountering an unrecognized string, continue to the next loop.
            uint32_t refLen, qryLen;
            if (!get_len(INFOSTRUCTTMP, VCFOPENCLASS, refLen, qryLen)) {
                continue;
            }

            uint32_t refEnd;
            refEnd = INFOSTRUCTTMP.POS + refLen - 1;

            // Determine the variation type and save it to the hash table
            if (refLen >= 50 || qryLen >= 50) {  // Save sv termination and mutation information
                svInfoMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS].information = INFOSTRUCTTMP.line;
                svInfoMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS].refEnd = refEnd;
            }
        }
    }


    // Find the number of SNPS and indel near sv and split them
    // Loop over the sv hash table
    for (const auto& [chromosome, refStartInfoMap] : svInfoMap) {
        // snp+indel start and end list
        const vector<uint32_t>& snpIndelRefStartVec = snpIndelInfoMap[chromosome].refStartVec;
        const vector<uint32_t>& snpIndelRefEndVec = snpIndelInfoMap[chromosome].refEndVec;

        for (const auto& [refStart, Info] : refStartInfoMap) {
            // Binary search finds the number of SNPS and indel within 200bp of starting and ending
            // Find left and right index
            int64_t startLeftIdxTmp = search_Binary_right(snpIndelRefEndVec, refStart - length_);
            int64_t startRightIdxTmp = search_Binary_left(snpIndelRefEndVec, refStart);
            int64_t endLeftIdxTmp = search_Binary_right(snpIndelRefStartVec, Info.refEnd);
            int64_t endRightIdxTmp = search_Binary_left(snpIndelRefStartVec, Info.refEnd + length_);

            // Number of snp+indel at both ends of sv
            int64_t leftSnpIndelNum = startRightIdxTmp - startLeftIdxTmp;
            if (leftSnpIndelNum < 0) {  // none
                leftSnpIndelNum = 0;
            } else {  // There is one or more
                leftSnpIndelNum++;
            }
            
            int64_t rightSnpIndelNum = endRightIdxTmp - endLeftIdxTmp;
            if (rightSnpIndelNum < 0) {  // none
                rightSnpIndelNum = 0;
            } else {  // There is one or more
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

            // save result
            if (SVs0Stream.tellp() >= CACHE_SIZE || 
                SVs1Stream.tellp() >= CACHE_SIZE || 
                SVs2Stream.tellp() >= CACHE_SIZE || 
                SVs4Stream.tellp() >= CACHE_SIZE || 
                SVs6Stream.tellp() >= CACHE_SIZE || 
                SVs8Stream.tellp() >= CACHE_SIZE || 
                SVsMoreStream.tellp() >= CACHE_SIZE
            ) {  // Determine the variation type and save it to the string
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

                // Clear stringstream
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
    
    // last save
    if (SVs0Stream.tellp() > 0 || 
        SVs1Stream.tellp() > 0 || 
        SVs2Stream.tellp() > 0 || 
        SVs4Stream.tellp() > 0 || 
        SVs6Stream.tellp() > 0 || 
        SVs8Stream.tellp() > 0 || 
        SVsMoreStream.tellp() > 0
    ) {  // Determine the variation type and save it to the string
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

        // Clear stringstream
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