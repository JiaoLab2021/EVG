#pragma once

#include <cstdint>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <future>
#include <memory>
#include <string>
#include <string_view>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iterator>
#include <unordered_map>

#include "hts_file.hpp"
#include "logger.hpp"
#include "task_pool.hpp"
#include "strip_split_join.hpp"

using namespace std;

// One VCF record split into its columns.
struct VCFINFOSTRUCT
{
    public:
    string line;
    vector<string> lineVec;
    vector<string> ALTVec;

    string CHROM;
    uint32_t POS;
    string ID;
    string REF;
    string ALT;
    double QUAL;
    string FILTER;
    string INFO;
    string FORMAT;

    uint32_t LEN;  // REF length
    uint32_t END;  // POS + LEN - 1

    VCFINFOSTRUCT() : POS(0), QUAL(0.0), LEN(0), END(0) {}

    // Reset every field, keeping string capacity for reuse.
    void clear() {
        line.clear();
        lineVec.clear();
        ALTVec.clear();
        CHROM.clear();
        POS = 0;
        ID.clear();
        REF.clear();
        ALT.clear();
        QUAL = 0.0;
        FILTER.clear();
        INFO.clear();
        FORMAT.clear();
        LEN = 0;
        END = 0;
    }

    // Exchange every field with other.
    void swap(VCFINFOSTRUCT& other) {
        using std::swap;
        swap(line, other.line);
        swap(lineVec, other.lineVec);
        swap(ALTVec, other.ALTVec);
        swap(CHROM, other.CHROM);
        swap(POS, other.POS);
        swap(ID, other.ID);
        swap(REF, other.REF);
        swap(ALT, other.ALT);
        swap(QUAL, other.QUAL);
        swap(FILTER, other.FILTER);
        swap(INFO, other.INFO);
        swap(FORMAT, other.FORMAT);
        swap(LEN, other.LEN);
        swap(END, other.END);
    }

    // Copy every field from other.
    VCFINFOSTRUCT& operator=(const VCFINFOSTRUCT& other) {
        if (this != &other) {
            line = other.line;
            lineVec = other.lineVec;
            ALTVec = other.ALTVec;
            CHROM = other.CHROM;
            POS = other.POS;
            ID = other.ID;
            REF = other.REF;
            ALT = other.ALT;
            QUAL = other.QUAL;
            FILTER = other.FILTER;
            INFO = other.INFO;
            FORMAT = other.FORMAT;
            LEN = other.LEN;
            END = other.END;
        }
        return *this;
    }

    // Move every field from other.
    VCFINFOSTRUCT(VCFINFOSTRUCT&& other) noexcept
        : line(std::move(other.line)),
        lineVec(std::move(other.lineVec)),
        ALTVec(std::move(other.ALTVec)),
        CHROM(std::move(other.CHROM)),
        POS(other.POS),
        ID(std::move(other.ID)),
        REF(std::move(other.REF)),
        ALT(std::move(other.ALT)),
        QUAL(other.QUAL),
        FILTER(std::move(other.FILTER)),
        INFO(std::move(other.INFO)),
        FORMAT(std::move(other.FORMAT)),
        LEN(other.LEN),
        END(other.END) {
    }

    // Move every field from other.
    VCFINFOSTRUCT& operator=(VCFINFOSTRUCT&& other) noexcept {
        if (this != &other) {
            line = std::move(other.line);
            lineVec = std::move(other.lineVec);
            ALTVec = std::move(other.ALTVec);
            CHROM = std::move(other.CHROM);
            POS = other.POS;
            ID = std::move(other.ID);
            REF = std::move(other.REF);
            ALT = std::move(other.ALT);
            QUAL = other.QUAL;
            FILTER = std::move(other.FILTER);
            INFO = std::move(other.INFO);
            FORMAT = std::move(other.FORMAT);
            LEN = other.LEN;
            END = other.END;
        }
        return *this;
    }
};

// A plain, gzip or BGZF VCF file read in chunks of whole records, the record parser and the static record helpers (worker-safe).
// Records are the lines read the way gzgets did: the final '\n' dropped ('\r' kept), a blank line joined to the next line.
class VCFFile
{
public:
    // Open the file; with threads > 1, BGZF input is inflated on the shared htslib pool. Throw if it cannot be opened.
    VCFFile(
        const string & vcfFileName,
        uint32_t threads
    );

    bool read_chunk(
        string & chunk
    );

    static vector<string_view> split_lines(
        const string & chunk
    );

    void parse_record(
        string_view line,
        VCFINFOSTRUCT & INFOSTRUCTTMP
    ) const;

    static string get_TYPE(
        const uint32_t & refLen,
        const vector<string> & ALTVec
    );

    static tuple<unordered_map<int, vector<string> >, uint32_t> get_gt(
        const vector<string> & lineVec
    );

    static vector<string> gt_split(string & gtTxt);

    static string sample_gt(
        const vector<string> & lineVec,
        size_t column,
        uint32_t gtIndex
    );

    static tuple<double, double> calculate(
        const vector<string> & lineVec,
        uint32_t misSampleNum
    );

protected:
    string vcfFileName_;
    uint32_t threads_;

private:
    InputFile file_;
    string pending_;  // text read after the last chunk boundary
    bool eof_;
};

// The result type of a reader without a per-record worker.
struct NoWork {};

// One record from a VCFReader: its columns, the per-record worker's result (value), and what the worker printed and threw.
// Records are reused for later chunks; the worker of each data record starts from a default Result.
template <class Result>
struct ProcessedRecord : Outcome<Result> {
    VCFINFOSTRUCT info;

    // Exchange the record with other, keeping string storage in circulation.
    void swap(ProcessedRecord& other) {
        using std::swap;
        info.swap(other.info);
        swap(this->value, other.value);
        swap(this->error, other.error);
        this->diagnostics.swap(other.diagnostics);
    }
};

// Read a VCF file record by record. Chunks of records are parsed ahead, 2 * threads at a time on the shared pool (inline with
// threads <= 1); with a worker, work(info, result) also runs there for every data record (no '#'), its output and exception kept
// in the record for the caller to replay with get(). A bad line or a read error is thrown when the caller reaches its position.
// Call read() from the main thread only, never from a pool task.
template <class Result = NoWork>
class VCFReader : public VCFFile
{
public:
    using Work = function<void(VCFINFOSTRUCT&, Result&)>;

    // Open the file; work, if given, runs on every data record.
    VCFReader(
        const string & vcfFileName,
        uint32_t threads = 1,
        Work work = nullptr
    ) : VCFFile(vcfFileName, threads), work_(std::move(work)), nextRecord_(0) {}

    // Wait for the chunks still being parsed.
    ~VCFReader() {
        for (auto& chunk : parsing_) {
            if (chunk.valid()) chunk.wait();
        }
    }

    // Swap the next record into record (the caller's previous record is reused later); false (record cleared) at the end of the file.
    bool read(ProcessedRecord<Result> & record) {
        if (!next_record()) {
            record.info.clear();
            return false;
        }
        record.swap(records_[nextRecord_++]);
        return true;
    }

    // Swap the next record's columns into INFOSTRUCTTMP; false (cleared) at the end of the file.
    bool read(VCFINFOSTRUCT & INFOSTRUCTTMP) {
        if (!next_record()) {
            INFOSTRUCTTMP.clear();
            return false;
        }
        INFOSTRUCTTMP.swap(records_[nextRecord_++].info);
        return true;
    }

private:
    // Records of one chunk, parsed up to the first bad line, whose exception is kept.
    struct ParsedChunk {
        vector<ProcessedRecord<Result> > records;
        exception_ptr error;
    };

    Work work_;
    deque<future<ParsedChunk> > parsing_;          // chunks being parsed, oldest first
    vector<ProcessedRecord<Result> > records_;     // parsed records not yet returned
    vector<ProcessedRecord<Result> > spare_;       // returned records, refilled by the next chunk so their strings keep their buffers
    size_t nextRecord_;
    exception_ptr badRecord_;                       // the bad line after records_
    exception_ptr readError_;                       // a read error, raised after the chunks read before it

    // Make records_[nextRecord_] the next record, keeping chunks parsing ahead; rethrow a bad line or read error at its position; false at the end.
    bool next_record() {
        while (nextRecord_ == records_.size()) {
            if (badRecord_) rethrow_exception(badRecord_);

            // ------ Keep chunks parsing ahead (the future's slot exists before its task does) ------
            const size_t maxParsing = threads_ > 1 ? 2 * static_cast<size_t>(threads_) : 1;
            while (!readError_ && parsing_.size() < maxParsing) {
                string chunk;
                try {
                    if (!read_chunk(chunk)) break;
                } catch (...) {
                    readError_ = current_exception();
                    break;
                }
                parsing_.emplace_back();
                try {
                    parsing_.back() = start_task(threads_, [this, chunk = std::move(chunk), records = std::move(spare_)]() mutable {
                        return parse_chunk(chunk, std::move(records));
                    });
                } catch (...) {
                    parsing_.pop_back();
                    throw;
                }
                spare_.clear();
            }

            // ------ Take the oldest parsed chunk ------
            if (parsing_.empty()) {
                if (readError_) rethrow_exception(readError_);
                return false;
            }
            ParsedChunk parsed = parsing_.front().get();
            parsing_.pop_front();
            spare_ = std::move(records_);
            records_ = std::move(parsed.records);
            nextRecord_ = 0;
            badRecord_ = parsed.error;
        }
        return true;
    }

    // Parse a chunk into recycled records, run the worker on data records with their output sliced per record; stop at a bad line.
    ParsedChunk parse_chunk(
        const string & chunk,
        vector<ProcessedRecord<Result> > records
    ) const {
        ParsedChunk parsed;
        const vector<string_view> lines = split_lines(chunk);
        parsed.records = std::move(records);
        parsed.records.resize(lines.size());

        unique_ptr<DiagnosticCapture> capture;
        vector<size_t> ends;
        if (work_) {
            capture = make_unique<DiagnosticCapture>();
            ends.reserve(lines.size());
        }

        for (size_t i = 0; i < lines.size(); ++i) {
            ProcessedRecord<Result>& record = parsed.records[i];
            if (record.info.line.capacity() > 2 * lines[i].size() + 4096) {
                ProcessedRecord<Result>().swap(record);  // free a record that once held a much longer line
            }
            try {
                parse_record(lines[i], record.info);
            } catch (...) {
                parsed.records.resize(i);
                parsed.error = current_exception();
                break;
            }
            record.error = nullptr;
            if (work_) {
                record.value = Result();  // records are recycled: the worker starts from a fresh result
                if (record.info.line.find('#') == string::npos) {
                    try {
                        work_(record.info, record.value);
                    } catch (...) {
                        record.error = current_exception();
                    }
                }
                ends.push_back(capture->size());
            }
        }

        if (work_) {
            const string& text = capture->text();
            size_t begin = 0;
            for (size_t i = 0; i < parsed.records.size(); ++i) {
                parsed.records[i].diagnostics.assign(text, begin, ends[i] - begin);
                begin = ends[i];
            }
        }
        return parsed;
    }
};

// The reader without a per-record worker.
using VCFOPEN = VCFReader<>;

// Throw if the records of a VCF file are not sorted by position within each chromosome.
void check_vcf_sort(const string & inputVcf);

// Return the REF length and the ALT lengths of a record: every ALT allele ("all") or the alleles in gtVec ("hap").
tuple<uint32_t, vector<uint32_t> > get_hap_len(
    const string & svType,
    const string & refSeq,
    const string & qrySeqs,
    const vector<int> & gtVec,
    const string & lenType
);
