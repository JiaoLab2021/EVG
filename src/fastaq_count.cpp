#include "../include/fastaq_count.hpp"

#include <exception>
#include <future>
#include <stdexcept>
#include <string>
#include <vector>

#include <htslib/kseq.h>

#include "../include/hts_file.hpp"
#include "../include/save.hpp"
#include "../include/task_pool.hpp"

KSEQ_INIT(InputFile*, read_input)

namespace {

constexpr uint64_t batch_bases = 1 << 20;

// Record lengths read from one file and what ended the batch: the kseq_read status (>= 0: more follow, -1: end of
// file, < -1: damaged record) or a read error, raised only after the lengths read before it.
struct ReadBatch {
    std::vector<uint64_t> lengths;
    int status = 0;
    std::exception_ptr error;
};

// Read records until they hold batch_bases bases, the file ends or a read error stops it.
ReadBatch read_batch(kseq_t* ks) {
    ReadBatch batch;
    uint64_t bases = 0;
    try {
        while (bases < batch_bases && (batch.status = kseq_read(ks)) >= 0) {
            batch.lengths.push_back(ks->seq.l);
            bases += ks->seq.l;
        }
    } catch (...) {
        batch.error = std::current_exception();
        batch.status = -1;
    }
    return batch;
}

// Record lengths of one file, read one batch ahead on the shared pool so two files are read at the same time.
class BatchReader {
public:
    BatchReader(kseq_t* ks, uint32_t threads) : ks_(ks), threads_(threads) { next_ = read_ahead(); }

    ~BatchReader() {
        if (next_.valid()) next_.wait();
    }

    BatchReader(const BatchReader&) = delete;
    BatchReader& operator=(const BatchReader&) = delete;

    // Give the next record length; false at end of file or on a damaged record, see status(); throw a read error
    // when its position is reached.
    bool next(uint64_t& length) {
        while (index_ == batch_.lengths.size()) {
            if (batch_.error) std::rethrow_exception(batch_.error);
            if (batch_.status < 0) return false;
            batch_ = next_.get();
            index_ = 0;
            if (batch_.status >= 0) next_ = read_ahead();
        }
        length = batch_.lengths[index_++];
        return true;
    }

    // The kseq_read status that ended the file once next() returned false.
    int status() const { return batch_.status; }

private:
    // Start reading the next batch.
    std::future<ReadBatch> read_ahead() {
        return start_task(threads_, [ks = ks_] { return read_batch(ks); });
    }

    kseq_t* ks_;
    uint32_t threads_;
    ReadBatch batch_;
    size_t index_ = 0;
    std::future<ReadBatch> next_;
};

}

// Store the input and output file names and the thread count.
Count::Count(const std::string& inputFileName1, const std::string& inputFileName2, const std::string& outputFileName,
             uint32_t threads)
    : inputFileName1_(inputFileName1), inputFileName2_(inputFileName2), outputFileName_(outputFileName),
      threads_(threads) {}

// Count reads and bases; paired files are read concurrently and zipped in step until either one ends.
void Count::count_reads() {
    if (inputFileName1_.length() > 0 && inputFileName2_.length() > 0) {
        // ------ Paired-end reads ------
        InputFile file1(inputFileName1_, threads_);
        InputFile file2(inputFileName2_, threads_);
        kseq_t* ks1 = kseq_init(&file1);
        kseq_t* ks2 = kseq_init(&file2);
        int status1 = 0;
        int status2 = 0;
        {
            BatchReader reader1(ks1, threads_);
            BatchReader reader2(ks2, threads_);
            uint64_t length1 = 0;
            uint64_t length2 = 0;
            while (true) {
                if (!reader1.next(length1)) {
                    status1 = reader1.status();
                    break;
                }
                if (!reader2.next(length2)) {
                    status2 = reader2.status();
                    break;
                }
                countOut_.readBase += length1;
                countOut_.readBase += length2;
                countOut_.readNum += 2;
            }
        }

        kseq_destroy(ks1);
        kseq_destroy(ks2);
        if (status1 < -1 || status2 < -1) {
            throw std::runtime_error("'" + (status1 < -1 ? inputFileName1_ : inputFileName2_)
                                     + "' -> truncated or damaged FASTA/FASTQ file.");
        }
    } else if (inputFileName1_.length() > 0 && inputFileName2_.length() == 0) {
        // ------ Single file ------
        InputFile file(inputFileName1_, threads_);
        kseq_t* ks = kseq_init(&file);
        int status = 0;
        while ((status = kseq_read(ks)) >= 0) {
            countOut_.readBase += ks->seq.l;
            countOut_.readNum++;
        }

        kseq_destroy(ks);
        if (status < -1) {
            throw std::runtime_error("'" + inputFileName1_ + "' -> truncated or damaged FASTA/FASTQ file.");
        }
    }
}

// Write readBase, readNum and readLen (integer mean, 0 for empty input) to the output file or stdout.
void Count::save_result() {
    const uint64_t readLen = countOut_.readNum > 0 ? countOut_.readBase / countOut_.readNum : 0;
    std::string outTxt = "readBase:" + std::to_string(countOut_.readBase) + "\n" +
                         "readNum:" + std::to_string(countOut_.readNum) + "\n" +
                         "readLen:" + std::to_string(readLen) + "\n";

    SAVE saveClass(outputFileName_, threads_);
    saveClass.save(outTxt);
}
