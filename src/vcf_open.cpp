#include "../include/vcf_open.hpp"

using namespace std;


/**
 * @brief initialize
 * 
 * @param vcfFileName   the output of vcf file
 * 
**/
VCFOPEN::VCFOPEN(
    const string & vcfFileName
) {
    vcfFileName_ = vcfFileName;

    // Input file stream
    gzfpI = gzopen(vcfFileName_.c_str(), "rb");
    if(!gzfpI) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << vcfFileName_ << "': No such file or directory." << endl;
        exit(1);
    }

    // buffer size
    bufferSize_ = 1 * 1024 * 1024;  // 1 Mb
    line_ = new char[bufferSize_];  // Read only bufferSize_ bytes of data at a time
}

/**
 * @brief Close file
 * 
 * @return 0
**/
VCFOPEN::~VCFOPEN() {
    // Close file
    gzclose(gzfpI);

    // Release dynamically allocated memory
    delete[] line_;
}


/**
 * @brief traverse files
 * 
 * @param INFOSTRUCTTMP store the contents of the line
 * 
 * @return bool
**/
bool VCFOPEN::read(
    VCFINFOSTRUCT & INFOSTRUCTTMP
) {
    INFOSTRUCTTMP.clear();

    string info = ""; // temporary string

    if(gzgets(gzfpI, line_, bufferSize_)) {
        info += line_;

        // empty line_
        memset(line_, 0, min(strlen(line_), bufferSize_));

        // If there is no newline character, it means that the line is not over, continue to read
        while ((info.find("\n") == string::npos || info == "\n") && gzgets(gzfpI, line_, bufferSize_)) {
            info += line_;

            // empty line_
            memset(line_, 0, min(strlen(line_), bufferSize_));
        }
        
        // remove line breaks
        if (!info.empty()) {
            info = strip(info, '\n');
        } else {  // false
            return false;
        }
        
        // assignment
        INFOSTRUCTTMP.line = info;
        INFOSTRUCTTMP.lineVec = split(info, "\t");  // segmentation

        // comment lines
        if (info.find("#") != string::npos) {
            return true;
        }
        

        // non-sample lines
        if (INFOSTRUCTTMP.lineVec.size() < 9) { // Check the file first, if it is wrong, it will jump out of the code
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: '"
                << vcfFileName_ 
                << "' -> The number of vcf columns is less than 9.\n" 
                << INFOSTRUCTTMP.line << endl;
            exit(1);
        }

        INFOSTRUCTTMP.CHROM = INFOSTRUCTTMP.lineVec[0];  // CHROM
        INFOSTRUCTTMP.POS = stoul(INFOSTRUCTTMP.lineVec[1]);  // POS
        INFOSTRUCTTMP.ID = INFOSTRUCTTMP.lineVec[2];  // ID
        INFOSTRUCTTMP.REF = INFOSTRUCTTMP.lineVec[3];  // REF
        INFOSTRUCTTMP.ALT = INFOSTRUCTTMP.lineVec[4];  // ALT
        INFOSTRUCTTMP.ALTVec = split(INFOSTRUCTTMP.ALT, ",");  // 'ALT' 'vec'

        if (isdigit(INFOSTRUCTTMP.lineVec[5][0])) {  // Assign QUAL if it is a number, otherwise 0.0
            INFOSTRUCTTMP.QUAL = stod(INFOSTRUCTTMP.lineVec[5]);  // QUAL
        } else {
            INFOSTRUCTTMP.QUAL = 0.0;
        }

        INFOSTRUCTTMP.FILTER = INFOSTRUCTTMP.lineVec[6];  // FILTER
        INFOSTRUCTTMP.INFO = INFOSTRUCTTMP.lineVec[7];  // [7]
        INFOSTRUCTTMP.FORMAT = INFOSTRUCTTMP.lineVec[8];  // [8]

        INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.REF.size();  // REF length
        INFOSTRUCTTMP.END = INFOSTRUCTTMP.POS + INFOSTRUCTTMP.LEN - 1;  // end position

        return true;
    } else {  // Return false directly after traversing
        return false;
    }
}


/**
 * Get the type of variant
 * 
 * @param refLen  the length of REF variant
 * @param ALTVec  vector<qrySeq>
 * 
 * @return string   TYPE: SNP, InDel, Deletion, Insertion, Inversion, Duplication, Other
**/
string VCFOPEN::get_TYPE(
    const uint32_t & refLen,
    const vector<string> & ALTVec
) {
    uint32_t qryLen = 0;
    for (const auto& it1 : ALTVec) {
        uint32_t qryLenTmp = it1.size();
        qryLen = max(qryLen, qryLenTmp);  // Find the longest qrySeq
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
    
    return "Other";
}


/**
 * Get a list of locus genotypes.
 *
 * @param lineTmpVec  lineTmpVec
 * 
 * 
 * @return GTVecMap, misSampleNum   map<idx, vector<GTString> >
**/
tuple<unordered_map<int, vector<string> >, uint32_t> VCFOPEN::get_gt(
    const vector<string> & lineTmpVec
) {
    unordered_map<int, vector<string> > GTVecMap;  // map of all Line types, map<idx, vector<GTString> >
    uint32_t misSampleNum = 0;  // number of missing samples

    // FORMAT column
    vector<string> formatVec = split(lineTmpVec[8], ":");  // FORMAT character split

    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // Gets the index location of GT

    if (gtIndex == formatVec.size()) {  // Check whether index exists. If no index exists, the genotypes are all 0.
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << lineTmpVec[0] << ":" << lineTmpVec[1] << endl;
        for (size_t i = 9; i < lineTmpVec.size(); i++) {  // Start the loop in the ninth column
            GTVecMap[i - 9] = {"0", "0"};  // Assign
        }
    } else {  // If it exists, save it
        string GTString;  // Store the genotype field

        for (size_t i = 9; i < lineTmpVec.size(); i++) {  // Start the loop in the ninth column
            GTString = split(lineTmpVec[i], ":")[gtIndex];  // gt field
            vector<string> GTVec = gt_split(GTString);
            GTVecMap[i - 9] = GTVec;  // Assign

            // Determine whether it is a missing sample
            if (GTString == "." || GTString == "./." || GTString == ".|." || GTString == "././." || GTString == ".|.|." || GTString == "./././." || GTString == ".|.|.|.") {
                misSampleNum++;
            }
        }
    }

    return tuple(GTVecMap, misSampleNum);
}


/**
 * split gt
 *
 * @param gtTxt
 * 
 * 
 * @return vector<gt>
**/
vector<string> VCFOPEN::gt_split(
    string & gtTxt
) {
    // Temporary gt list
    vector<string> gtVecTmp;

    // Replace all '.' with '0'
    std::replace(gtTxt.begin(), gtTxt.end(), '.', '0');
    
    if (gtTxt.find("/") != string::npos) {
        gtVecTmp = split(gtTxt, "/");
    } else if (gtTxt.find("|") != string::npos) {
        gtVecTmp = split(gtTxt, "|");
    } else {
        try {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: sample has only one genotype, attempting to correct to diploid -> " << gtTxt << endl;
            stoul(gtTxt);
            gtVecTmp.push_back(gtTxt);
        } catch (const std::invalid_argument&) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: GT is not separated by '/' or '|' -> " << gtTxt << endl;
            exit(1);
        } catch (const std::out_of_range&) {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: GT is not separated by '/' or '|' -> " << gtTxt << endl;
            exit(1);
        }
    }

    return gtVecTmp;
}


/**
 * Calculate MAF and MISSRATE.
 *
 * @param GTVecMap      map<int, vector<string> >,  map<idx, vector<GTString> >
 * @param misSampleNum  number of missing samples
 * 
 * 
 * @return tuple<double, double>   tuple<MAF, MISSRATE>
**/
tuple<double, double> VCFOPEN::calculate(
    const unordered_map<int, vector<string> > & GTVecMap, 
    uint32_t misSampleNum
) {
    double MAFTMP;  // Minimum allele frequency
    double MISSRATETMP;  // Miss rate

    uint32_t sampleNum = GTVecMap.size();

    /* ******************************** Calculate MAF ******************************** */ 
    map<string, uint32_t> GTFreMapTmp;  // Genotype frequency map<gt, fre>
    uint16_t ploidy = 1;  // Record ploidy

    for (const auto& it1 : GTVecMap) {  // map<idx, vector<GTString> >
        for (const auto& it2 : it1.second) {  // vector<GTString>
            GTFreMapTmp[it2]++;  // Add 1 to the frequency corresponding to the genotype
        }
        // If it is 1 then calculate again
        if (ploidy == 1) {
            ploidy = it1.second.size();  // Record ploidy
        }
    }

    map<uint32_t, string> GTFreMapTmpTmp;  // map conversion, map<fre, gt>
    for (const auto& it1 : GTFreMapTmp) {  // map<gt, fre>
        GTFreMapTmpTmp[it1.second] = it1.first;
    }

    if (GTFreMapTmpTmp.size() > 1) {
        // The second-to-last is the suballele
        auto iter1 = GTFreMapTmpTmp.end();
        iter1--; iter1--;
        // frequence/(ploidy*N)
        MAFTMP = iter1->first/(double)(ploidy*sampleNum);
    }
    
    /* ******************************** Calculate MISSRATE ******************************** */ 
    // 1- number/allNumber
    MISSRATETMP = 1 - (misSampleNum/(double)sampleNum);

    return make_tuple(MAFTMP, MISSRATETMP);
}