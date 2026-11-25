#include "../include/vcf_breakpoint.hpp"

#include "../include/logger.hpp"

using namespace std;

// Store the input VCF, output prefix, error size and thread count.
VCFBreakpoint::VCFBreakpoint(
    const string & vcfFileName,
    const string & prefix,
    const int32_t & breakpointErrorSize,
    uint32_t threads
) : vcfFileName_(vcfFileName), prefix_(prefix), breakpointErrorSize_(breakpointErrorSize), threads_(threads) {}

// Subtract the error size from POS (minimum 1), update END in INFO and write <prefix>.vcf.gz; workers shift the records (shift_record).
void VCFBreakpoint::shift_breakpoints()
{
    VCFReader<string> VCFOPENCLASS(vcfFileName_, threads_, [this](VCFINFOSTRUCT& INFOSTRUCTTMP, string& outLine) {
        shift_record(INFOSTRUCTTMP, outLine);
    });

    string outTxt;

    SAVE outFile(prefix_ + ".vcf.gz", threads_);

    ProcessedRecord<string> record;
    while (VCFOPENCLASS.read(record)) {
        if (record.info.line.find("#") != string::npos) {
            string headLine = record.info.line + "\n";
            outFile.save(headLine);
        } else {
            outTxt += record.get();

            if (outTxt.size() > 10 * 1024 * 1024) {
                outFile.save(outTxt);

                outTxt.clear();
            }
        }
    }

    outFile.save(outTxt);

    string().swap(outTxt);
}

// Shift one data record and set outLine to the record as written (joined, '\n').
void VCFBreakpoint::shift_record(VCFINFOSTRUCT & INFOSTRUCTTMP, string & outLine) const
{
    static const std::regex reg("END=(\\d+)");

    INFOSTRUCTTMP.lineVec[1] = (static_cast<int32_t>(INFOSTRUCTTMP.POS) > breakpointErrorSize_) ? to_string(static_cast<int32_t>(INFOSTRUCTTMP.POS) - breakpointErrorSize_) : "1";

    uint32_t refEnd = stoul(INFOSTRUCTTMP.lineVec[1]) + INFOSTRUCTTMP.LEN - 1;
    INFOSTRUCTTMP.lineVec[7] = regex_replace(string(INFOSTRUCTTMP.lineVec[7]), reg, string("END=" + to_string(refEnd)));

    outLine = join(INFOSTRUCTTMP.lineVec, "\t");
    outLine += '\n';
}
