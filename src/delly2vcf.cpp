#include "../include/delly2vcf.hpp"

#include <stdexcept>

#include <htslib/kseq.h>

#include "../include/hts_file.hpp"
#include "../include/logger.hpp"

using namespace std;

KSEQ_INIT(InputFile*, read_input)

// Store the reference, DELLY VCF and output names and the thread count.
Delly2VCF::Delly2VCF(
    const string& refFileName,
    const string& dellyFilename,
    const string& outputFilename,
    uint32_t threads
) : refFileName_(refFileName), dellyFilename_(dellyFilename), outputFilename_(outputFilename), threads_(threads) {}

// Load the reference sequences into memory.
void Delly2VCF::build_reference_index()
{
    log_stream() << "Building refFileName_ index: " << refFileName_ << endl;

    InputFile reference(refFileName_, threads_);
    kseq_t *ks = kseq_init(&reference);

    int status;
    while( (status = kseq_read(ks)) >= 0 ) {
        FastaMap_.emplace(ks->name.s, ks->seq.s);
    }
    if (status < -1) {
        throw runtime_error("'" + refFileName_ + "' -> truncated or damaged FASTA/FASTQ file.");
    }
    kseq_destroy(ks);
}

// Write the records converted by the workers (convert_record) in input order.
void Delly2VCF::convert_vcf()
{
    SAVE SAVEClass(outputFilename_, threads_);

    stringstream outStream;
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10;
    outStream.str().reserve(CACHE_SIZE);

    VCFReader<string> VCFOPENCLASS(dellyFilename_, threads_, [this](VCFINFOSTRUCT& INFOSTRUCTTMP, string& outLine) {
        convert_record(INFOSTRUCTTMP, outLine);
    });

    ProcessedRecord<string> record;
    while (VCFOPENCLASS.read(record)) {
        if (record.info.line.empty()) {
            continue;
        }

        if (record.info.line.find("#") != string::npos) {
            outStream << record.info.line << endl;
            continue;
        }

        const string& outLine = record.get();
        if (outLine.empty()) {
            continue;
        }

        outStream << outLine;

        if (outStream.tellp() >= CACHE_SIZE) {
            string outTxt = outStream.str();
            SAVEClass.save(outTxt);
            outStream.str(string());
            outStream.clear();
        }
    }

    if (outStream.tellp() > 0) {
        string outTxt = outStream.str();
        SAVEClass.save(outTxt);
        outStream.str(string());
        outStream.clear();
    }
}

// Set outLine to a PASS record with END, at most 100 kb, REF from the reference and ALT built for <INV>/<DEL>/<DUP>; "" if skipped.
void Delly2VCF::convert_record(VCFINFOSTRUCT& INFOSTRUCTTMP, string& outLine) const
{
    static const std::regex end_reg("END=(\\d+)");
    static const std::regex reg1("\\[");
    static const std::regex reg2("\\]");

    if (INFOSTRUCTTMP.FILTER != "PASS") {
        return;
    }

    INFOSTRUCTTMP.ALT = regex_replace(string(INFOSTRUCTTMP.ALT), reg1, string(""));
    INFOSTRUCTTMP.ALT = regex_replace(string(INFOSTRUCTTMP.ALT), reg2, string(""));

    smatch end_search;
    uint32_t refEnd = 0;
    if (regex_search(INFOSTRUCTTMP.INFO, end_search, end_reg)) {
        refEnd = stoul(end_search.str(1));
    } else {
        warning_stream() << "INFO column is missing END information, skipping this position: " << INFOSTRUCTTMP.line << endl;
        return;
    }

    uint32_t refLen = refEnd - INFOSTRUCTTMP.POS + 1;
    if (refLen > 0.1 * 1000 * 1000) {
        warning_stream() << "The variant is too long, skipping this position: " << INFOSTRUCTTMP.line << endl;
        return;
    }

    const auto chromosome = FastaMap_.find(INFOSTRUCTTMP.CHROM);
    const string emptySequence;
    string refSeq = (chromosome == FastaMap_.end() ? emptySequence : chromosome->second).substr(INFOSTRUCTTMP.POS - 1, refLen);
    string qrySeq = "";

    if (INFOSTRUCTTMP.ALT == "<INV>") {
        qrySeq = refSeq;
        reverse(qrySeq.begin(), qrySeq.end());
    } else if (INFOSTRUCTTMP.ALT == "<DEL>") {
        qrySeq = refSeq[0];
    } else if (INFOSTRUCTTMP.ALT == "<DUP>") {
        qrySeq = refSeq + refSeq;
    } else if (INFOSTRUCTTMP.ALT == "<INS>") {
        warning_stream() << "ALT == <INS>, skipping this position: " << INFOSTRUCTTMP.line << endl;
        return;
    } else if (INFOSTRUCTTMP.ALT.find(":") != string::npos) {
        warning_stream() << "Skipping translocation: " << INFOSTRUCTTMP.line << endl;
        return;
    } else {
        qrySeq = INFOSTRUCTTMP.ALT;
    }

    INFOSTRUCTTMP.lineVec[3] = refSeq;
    INFOSTRUCTTMP.lineVec[4] = qrySeq;

    outLine = join(INFOSTRUCTTMP.lineVec, "\t");
    outLine += '\n';
}
