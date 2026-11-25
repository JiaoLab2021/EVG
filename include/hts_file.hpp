#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <sys/types.h>

struct BGZF;

// (De)compress BGZF on the shared htslib pool (`threads` threads, writers queue 256 blocks); return false for plain,
// gzip and threads <= 1, which stay single-threaded.
bool attach_bgzf_pool(BGZF* fp, uint32_t threads);

// True for a regular file (not a pipe, FIFO or device), following symlinks.
bool is_regular_file(const std::string& name);

// A plain, gzip or BGZF input file read as bytes. Any thread count gives the same bytes and stops at the same
// damaged block: htslib's threaded reader can fail early, so a failed threaded read is retried single-threaded
// (only regular files are read on threads, as the retry reopens the file).
class InputFile
{
public:
    // Open the file; with threads > 1, BGZF in a regular file is inflated on the shared htslib pool. Throw if it cannot be opened.
    InputFile(const std::string& fileName, uint32_t threads);

    ~InputFile();

    InputFile(const InputFile&) = delete;
    InputFile& operator=(const InputFile&) = delete;

    // Read up to size bytes; return the count, 0 at the end; throw on damaged input or a missing BGZF EOF marker.
    ssize_t read(void* buffer, size_t size);

private:
    // Reopen single-threaded and skip the bytes already returned.
    void reopen_single_threaded();

    std::string fileName_;
    BGZF* fp_;
    bool threaded_;
    uint64_t returned_ = 0;  // bytes returned so far
};

// The kseq read function: KSEQ_INIT(InputFile*, read_input).
inline ssize_t read_input(InputFile* file, void* buffer, size_t size) { return file->read(buffer, size); }
