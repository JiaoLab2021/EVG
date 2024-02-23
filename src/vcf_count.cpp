// g++ vcf_count.cpp -o vcf_count -lz -O2
#include "../include/vcf_count.hpp"


using namespace std;


int main_count(int argc, char** argv)
{
    // Input file
    string vcfFileName;

    // Input parameter
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf", required_argument, 0, 'v'},

            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "v:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'v':
            vcfFileName = optarg;
            break;
        case 'h':
        case '?':
            help_count(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // Check whether the parameters are correct
    if (vcfFileName.empty())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Parameter error.\n";
        help_count(argv);
        return 1;
    }

    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";
    
    // init
    VCFCount VCFCountClass(vcfFileName);
    // count
    VCFCountClass.count();
    // count number and length
    VCFCountClass.count_num_len();
    // result
    VCFCountClass.get_result();
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// Help document
void help_count(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -v [options]" << endl
       << "calculate the number of alleles." << endl
       << endl
       << "required arguments:" << endl
       << "    -v, --vcf           FILE       vcf file to be converted" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}


/**
 * init
 *
 * @param vcfFileName       input VCF  file name
 * 
**/
VCFCount::VCFCount(
    const string & vcfFileName
) : vcfFileName_(vcfFileName) {
    // Variable length range
    lengthVec_.push_back("-999999/-10000");
    lengthVec_.push_back("-9999/-5000");
    lengthVec_.push_back("-4999/-2500");
    lengthVec_.push_back("-2499/-1000");
    lengthVec_.push_back("-999/-500");
    lengthVec_.push_back("-499/-400");
    lengthVec_.push_back("-399/-300");
    lengthVec_.push_back("-299/-200");
    lengthVec_.push_back("-199/-100");
    lengthVec_.push_back("-99/-50");
    lengthVec_.push_back("-49/-1");
    lengthVec_.push_back("0/0");
    lengthVec_.push_back("1/49");
    lengthVec_.push_back("50/99");
    lengthVec_.push_back("100/199");
    lengthVec_.push_back("200/299");
    lengthVec_.push_back("300/399");
    lengthVec_.push_back("400/499");
    lengthVec_.push_back("500/999");
    lengthVec_.push_back("1000/2499");
    lengthVec_.push_back("2500/4999");
    lengthVec_.push_back("5000/9999");
    lengthVec_.push_back("10000/999999");
}

/**
 * count
 *
 * @return void
**/
void VCFCount::count() {
    // Store vcf information
    VCFINFOSTRUCT INFOSTRUCTTMP;

    // initialize
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // If not, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty or comment line, skip
        if (INFOSTRUCTTMP.line.empty() || INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }
        
        // Traverse the qry sequence list
        for (const auto& qrySeq : INFOSTRUCTTMP.ALTVec) {
            // Record the amount of variation
            uint32_t qryLen = qrySeq.size();
            int32_t svLen = qryLen - INFOSTRUCTTMP.LEN;
            double lengthRatio = qryLen / float(INFOSTRUCTTMP.LEN);

            // Allelic Type
            bool mulAllelicBool = (INFOSTRUCTTMP.ALT.find(",") != string::npos) ? true : false;

            // SNPs
            if (svLen == 0 && INFOSTRUCTTMP.LEN == 1 && qryLen == 1) {
                snpNum_++;
                snpMuAllelicNum_ += mulAllelicBool;
                snpBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            // Indels
            } else if (svLen <= 49 && svLen >= -49 && INFOSTRUCTTMP.LEN <= 49 && qryLen <= 49) {
                indelNum_++;
                indelMuAllelicNum_ += mulAllelicBool;
                indelBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            // Inversion
            } else if (svLen >= -2 && svLen <= 2 && INFOSTRUCTTMP.LEN > 49 && qryLen > 49 ) {
                invNum_++;
                invMuAllelicNum_ += mulAllelicBool;
                invBiAllelicNum_ += !mulAllelicBool;
                InvLenVec_.push_back(qryLen);
            // Duplication
            } else if (lengthRatio >= 1.8 && lengthRatio <= 2.2 && INFOSTRUCTTMP.LEN > 49 && qryLen > 49) {
                dupNum_++;
                dupMuAllelicNum_ += mulAllelicBool;
                dupBiAllelicNum_ += !mulAllelicBool;
                DupLenVec_.push_back(svLen);
            // Deletion
            } else if (svLen < 0) {
                delNum_++;
                delMuAllelicNum_ += mulAllelicBool;
                delBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            // Insertion
            } else if (svLen > 0) {
                insNum_++;
                insMuAllelicNum_ += mulAllelicBool;
                insBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            // Other
            } else {
                otherNum_++;
                otherMuAllelicNum_ += mulAllelicBool;
                otherBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            }
        }
    }
}

/**
 * count number and length
 *
 * @return void
**/
void VCFCount::count_num_len() {
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        // number
        int64_t x1 = std::stoll(split(lengthVec_[i], "/")[0]);
        int64_t x2 = std::stoll(split(lengthVec_[i], "/")[1]);

        // sv
        vector<uint64_t>::size_type result = count_if(svLenVec_.begin(), svLenVec_.end(), f_mod_c(x1, x2));
        svNumVec_.push_back(result);
        // total length
        int64_t totalLen = 0;
        for (const auto& it : svLenVec_) {
            if (it >= x1 && it <= x2) {
                totalLen += it;
            }
        }
        svTotalLenVec_.push_back(abs(totalLen));

        // Inversion
        result = count_if(InvLenVec_.begin(), InvLenVec_.end(), f_mod_c(x1, x2));
        InvNumVec_.push_back(result);
        // total length
        totalLen = 0;
        for (const auto& it : InvLenVec_) {
            if (it >= x1 && it <= x2) {
                totalLen += it;
            }
        }
        InvTotalLenVec_.push_back(abs(totalLen));

        // Duplication
        result = count_if(DupLenVec_.begin(), DupLenVec_.end(), f_mod_c(x1, x2));
        DupNumVec_.push_back(result);
        // total length
        totalLen = 0;
        for (const auto& it : DupLenVec_) {
            if (it >= x1 && it <= x2) {
                totalLen += it;
            }
        }
        DupTotalLenVec_.push_back(abs(totalLen));
    }
}

/**
 * print
 *
 * @return void
**/
void VCFCount::get_result() {
    // number
    cout << "#[Type]: Variants Bi-allelic Multi-allelic\n";
    cout << "[SNP]: " << snpNum_ << " " << snpBiAllelicNum_ << " " << snpMuAllelicNum_ << endl
         << "[Indels]: " << indelNum_ << " " << indelBiAllelicNum_ << " " << indelMuAllelicNum_ << endl
         << "[Insertion]: " << insNum_ << " " << insBiAllelicNum_ << " " << insMuAllelicNum_ << endl
         << "[Deletion]: " << delNum_  << " " << delBiAllelicNum_ << " " << delMuAllelicNum_<< endl
         << "[Inversion]: " << invNum_ << " " << invBiAllelicNum_ << " " << invMuAllelicNum_ << endl
         << "[Duplication]: " << dupNum_ << " " << dupBiAllelicNum_ << " " << dupMuAllelicNum_ << endl
         << "[Other]: " << otherNum_ << " " << otherBiAllelicNum_ << " " << otherMuAllelicNum_ << endl << endl;

    // SV number
    cout << "##SNPs+Indels+Deletion+Insertion\n#[Range]: Number\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << svNumVec_[i] << endl;
    }
    // SV total length
    cout << "\n#[Range]: Length\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << svTotalLenVec_[i] << endl;
    }
    cout << endl << endl;

    // Inversion number
    cout << "##Inversion\n#[Range]: Number\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << InvNumVec_[i] << endl;
    }
    // Inversion total length
    cout << "\n#[Range]: Length\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << InvTotalLenVec_[i] << endl;
    }
    cout << endl << endl;

    // Duplication number
    cout << "##Duplication\n#[Range]: Number\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << DupNumVec_[i] << endl;
    }
    // Duplication total length
    cout << "\n#[Range]: Length\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << DupTotalLenVec_[i] << endl;
    }
    cout << endl << endl;
}