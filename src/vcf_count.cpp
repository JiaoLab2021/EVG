#include "../include/vcf_count.hpp"

#include <cstdlib>

#include "../include/logger.hpp"

using namespace std;

// Store the input VCF, the thread count and the length bins.
VCFCount::VCFCount(
    const string & vcfFileName,
    uint32_t threads
) : vcfFileName_(vcfFileName), threads_(threads) {
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

// Classify every ALT allele as SNP, indel, insertion, deletion, inversion, duplication or other.
void VCFCount::count_variants() {
    VCFINFOSTRUCT INFOSTRUCTTMP;

    VCFOPEN VCFOPENCLASS(vcfFileName_, threads_);

    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        if (INFOSTRUCTTMP.line.empty() || INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        for (const auto& qrySeq : INFOSTRUCTTMP.ALTVec) {
            uint32_t qryLen = qrySeq.size();
            int32_t svLen = qryLen - INFOSTRUCTTMP.LEN;
            double lengthRatio = qryLen / float(INFOSTRUCTTMP.LEN);

            bool mulAllelicBool = (INFOSTRUCTTMP.ALT.find(",") != string::npos) ? true : false;

            if (svLen == 0 && INFOSTRUCTTMP.LEN == 1 && qryLen == 1) {
                snpNum_++;
                snpMuAllelicNum_ += mulAllelicBool;
                snpBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            } else if (svLen <= 49 && svLen >= -49 && INFOSTRUCTTMP.LEN <= 49 && qryLen <= 49) {
                indelNum_++;
                indelMuAllelicNum_ += mulAllelicBool;
                indelBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            } else if (svLen >= -2 && svLen <= 2 && INFOSTRUCTTMP.LEN > 49 && qryLen > 49 ) {
                invNum_++;
                invMuAllelicNum_ += mulAllelicBool;
                invBiAllelicNum_ += !mulAllelicBool;
                InvLenVec_.push_back(qryLen);
            } else if (lengthRatio >= 1.8 && lengthRatio <= 2.2 && INFOSTRUCTTMP.LEN > 49 && qryLen > 49) {
                dupNum_++;
                dupMuAllelicNum_ += mulAllelicBool;
                dupBiAllelicNum_ += !mulAllelicBool;
                DupLenVec_.push_back(svLen);
            } else if (svLen < 0) {
                delNum_++;
                delMuAllelicNum_ += mulAllelicBool;
                delBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            } else if (svLen > 0) {
                insNum_++;
                insMuAllelicNum_ += mulAllelicBool;
                insBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            } else {
                otherNum_++;
                otherMuAllelicNum_ += mulAllelicBool;
                otherBiAllelicNum_ += !mulAllelicBool;
                svLenVec_.push_back(svLen);
            }
        }
    }
}

// Count the alleles and sum their lengths in each length bin.
void VCFCount::count_length_bins() {
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        int64_t x1 = std::stoll(split(lengthVec_[i], "/")[0]);
        int64_t x2 = std::stoll(split(lengthVec_[i], "/")[1]);

        vector<uint64_t>::size_type result = count_if(svLenVec_.begin(), svLenVec_.end(), f_mod_c(x1, x2));
        svNumVec_.push_back(result);
        int64_t totalLen = 0;
        for (const auto& it : svLenVec_) {
            if (it >= x1 && it <= x2) {
                totalLen += it;
            }
        }
        svTotalLenVec_.push_back(abs(totalLen));

        result = count_if(InvLenVec_.begin(), InvLenVec_.end(), f_mod_c(x1, x2));
        InvNumVec_.push_back(result);
        totalLen = 0;
        for (const auto& it : InvLenVec_) {
            if (it >= x1 && it <= x2) {
                totalLen += it;
            }
        }
        InvTotalLenVec_.push_back(abs(totalLen));

        result = count_if(DupLenVec_.begin(), DupLenVec_.end(), f_mod_c(x1, x2));
        DupNumVec_.push_back(result);
        totalLen = 0;
        for (const auto& it : DupLenVec_) {
            if (it >= x1 && it <= x2) {
                totalLen += it;
            }
        }
        DupTotalLenVec_.push_back(abs(totalLen));
    }
}

// Print the type and length tables to stdout.
void VCFCount::print_result() {
    cout << "#[Type]: Variants Bi-allelic Multi-allelic\n";
    cout << "[SNP]: " << snpNum_ << " " << snpBiAllelicNum_ << " " << snpMuAllelicNum_ << endl
         << "[Indels]: " << indelNum_ << " " << indelBiAllelicNum_ << " " << indelMuAllelicNum_ << endl
         << "[Insertion]: " << insNum_ << " " << insBiAllelicNum_ << " " << insMuAllelicNum_ << endl
         << "[Deletion]: " << delNum_  << " " << delBiAllelicNum_ << " " << delMuAllelicNum_<< endl
         << "[Inversion]: " << invNum_ << " " << invBiAllelicNum_ << " " << invMuAllelicNum_ << endl
         << "[Duplication]: " << dupNum_ << " " << dupBiAllelicNum_ << " " << dupMuAllelicNum_ << endl
         << "[Other]: " << otherNum_ << " " << otherBiAllelicNum_ << " " << otherMuAllelicNum_ << endl << endl;

    cout << "##SNPs+Indels+Deletion+Insertion\n#[Range]: Number\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << svNumVec_[i] << endl;
    }
    cout << "\n#[Range]: Length\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << svTotalLenVec_[i] << endl;
    }
    cout << endl << endl;

    cout << "##Inversion\n#[Range]: Number\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << InvNumVec_[i] << endl;
    }
    cout << "\n#[Range]: Length\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << InvTotalLenVec_[i] << endl;
    }
    cout << endl << endl;

    cout << "##Duplication\n#[Range]: Number\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << DupNumVec_[i] << endl;
    }
    cout << "\n#[Range]: Length\n";
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        cout << "[" << lengthVec_[i] << "]: " << DupTotalLenVec_[i] << endl;
    }
    cout << endl << endl;
}
