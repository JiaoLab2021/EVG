#include "../include/vcf_filter.hpp"

#include "../include/logger.hpp"

using namespace std;

// Store the input VCF, output name, thresholds and thread count.
VCFFilter::VCFFilter(
    const string & vcfFileName,
    const string & outputFileName,
    const double & MAF,
    const double & MISSRATE,
    uint32_t threads
) : vcfFileName_(vcfFileName), outputFileName_(outputFileName), MAF_(MAF), MISSRATE_(MISSRATE), threads_(threads) {}

// Keep header lines, non-SNP records and SNPs that pass both thresholds; workers decide per record (filter_record).
void VCFFilter::filter_vcf()
{
    VCFReader<bool> VCFOPENCLASS(vcfFileName_, threads_, [this](VCFINFOSTRUCT& INFOSTRUCTTMP, bool& keep) {
        keep = filter_record(INFOSTRUCTTMP);
    });

    SAVE SAVECLASS(outputFileName_, threads_);

    string outTxt = "";

    ProcessedRecord<bool> record;
    while (VCFOPENCLASS.read(record)) {
        const VCFINFOSTRUCT& INFOSTRUCTTMP = record.info;

        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            outTxt += INFOSTRUCTTMP.line + "\n";
            continue;
        }

        if (!record.get()) {
            continue;
        }
        outTxt += INFOSTRUCTTMP.line + "\n";

        if (outTxt.size() >= 10 * 1024 * 1024) {
            SAVECLASS.save(
                outTxt
            );

            outTxt.clear();
        }
    }

    if (outTxt.size() > 0) {
        SAVECLASS.save(
            outTxt
        );

        string().swap(outTxt);
    }
}

// Return whether a data record is written (on a worker): non-SNP records always, SNPs with 2 called samples passing both thresholds.
bool VCFFilter::filter_record(VCFINFOSTRUCT & INFOSTRUCTTMP) const
{
    INFOSTRUCTTMP.ID = VCFOPEN::get_TYPE(
        INFOSTRUCTTMP.LEN,
        INFOSTRUCTTMP.ALTVec
    );

    if (INFOSTRUCTTMP.ID != "SNP") {
        return true;
    }

    double MAFTMP;
    double MISSRATETMP;

    unordered_map<int, vector<string> > GTVecMapTmp;
    uint32_t misSampleNum;
    tie(GTVecMapTmp, misSampleNum) = VCFOPEN::get_gt(
        INFOSTRUCTTMP.lineVec
    );

    if (GTVecMapTmp.size() - misSampleNum <= 1) {
        return false;
    }

    tie(MAFTMP, MISSRATETMP) = VCFOPEN::calculate(
        INFOSTRUCTTMP.lineVec,
        misSampleNum
    );

    return MAFTMP >= MAF_ && MISSRATETMP <= MISSRATE_;
}
