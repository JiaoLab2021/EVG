#include "../include/hts_file.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include <sys/stat.h>

#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>

using namespace std;

// Attach the shared htslib pool, created on first use with `threads` threads, to a BGZF reader or writer.
bool attach_bgzf_pool(BGZF* fp, uint32_t threads)
{
    if (threads <= 1 || !fp->is_compressed || fp->is_gzip) {
        return false;
    }
    static hts_tpool* pool = hts_tpool_init(static_cast<int>(threads));
    if (pool == nullptr || bgzf_thread_pool(fp, pool, fp->is_write ? 256 : 0) < 0) {
        throw runtime_error("failed to start the htslib thread pool.");
    }
    return true;
}

// True for a regular file (not a pipe, FIFO or device), following symlinks.
bool is_regular_file(const string& name) {
    struct stat info;
    return stat(name.c_str(), &info) == 0 && S_ISREG(info.st_mode);
}

// Open the file and attach the htslib pool to BGZF in a regular file.
InputFile::InputFile(const string& fileName, uint32_t threads) : fileName_(fileName) {
    fp_ = bgzf_open(fileName_.c_str(), "r");
    if (fp_ == nullptr) {
        throw runtime_error("'" + fileName_ + "': No such file or directory.");
    }
    try {
        threaded_ = is_regular_file(fileName_) && attach_bgzf_pool(fp_, threads);
    } catch (...) {
        bgzf_close(fp_);
        throw;
    }
}

// Close the file.
InputFile::~InputFile()
{
    bgzf_close(fp_);
}

// Read the next bytes; a failed threaded read is retried single-threaded, so the error is the one -t 1 reports.
ssize_t InputFile::read(void* buffer, size_t size)
{
    ssize_t bytes = bgzf_read(fp_, buffer, size);
    if (bytes < 0 && threaded_) {
        reopen_single_threaded();
        bytes = bgzf_read(fp_, buffer, size);
    }
    if (bytes < 0) {
        throw runtime_error("'" + fileName_ + "' -> read error (damaged or truncated gzip file).");
    }
    if (bytes == 0 && fp_->no_eof_block) {
        throw runtime_error("'" + fileName_ + "' -> truncated BGZF file (EOF marker missing).");
    }
    returned_ += bytes;
    return bytes;
}

// Reopen the file without threads and read past the bytes already returned.
void InputFile::reopen_single_threaded()
{
    bgzf_close(fp_);
    threaded_ = false;
    fp_ = bgzf_open(fileName_.c_str(), "r");
    if (fp_ == nullptr) {
        throw runtime_error("'" + fileName_ + "': No such file or directory.");
    }
    vector<char> skipped(1 << 16);
    for (uint64_t left = returned_; left > 0; ) {
        const ssize_t bytes = bgzf_read(fp_, skipped.data(), min<uint64_t>(left, skipped.size()));
        if (bytes <= 0) {
            throw runtime_error("'" + fileName_ + "' -> read error (damaged or truncated gzip file).");
        }
        left -= bytes;
    }
}
