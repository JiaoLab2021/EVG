#include "../include/save.hpp"

#include <cstdio>
#include <exception>
#include <stdexcept>

#include <htslib/bgzf.h>

#include "../include/hts_file.hpp"

using namespace std;

namespace {

// Remove a regular file (never a device such as /dev/null).
void remove_regular_file(const string & name)
{
    if (is_regular_file(name)) {
        std::remove(name.c_str());
    }
}

// The write error of an output file, or of stdout for an empty name.
runtime_error write_error(const string & name)
{
    return runtime_error((name.empty() ? string("stdout") : "'" + name + "'") + ": write error.");
}

}

// Open the output: BGZF compressed on `threads` htslib threads when the name contains .gz or .GZ, a plain file for other names, stdout for an empty name; throw if it cannot be created.
SAVE::SAVE(
    const string & outputFileName,
    uint32_t threads
) : outputFileName_(outputFileName), bgzfO_(nullptr), uncaughtAtOpen_(uncaught_exceptions())
{
    const string openError = "'" + outputFileName_ + "': No such file or directory or possibly reached the maximum open file limit. You can set 'ulimit -n' to a larger value to continue.";

    if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos) {
        bgzfO_ = bgzf_open(outputFileName_.c_str(), "w5");
        if (bgzfO_ == nullptr) {
            throw runtime_error(openError);
        }
        try {
            attach_bgzf_pool(bgzfO_, threads);
        } catch (...) {
            bgzf_close(bgzfO_);
            remove_regular_file(outputFileName_);
            throw;
        }
    } else if (outputFileName_.size() > 0) {
        fpO.open(outputFileName_, ios::out);
        if (!fpO) {
            throw runtime_error(openError);
        }
    }
}

// Close the output (flush stdout); if that fails (e.g. disk full) remove the file and throw; while an exception propagates, close quietly and remove the file.
SAVE::~SAVE() noexcept(false)
{
    const bool unwinding = uncaught_exceptions() > uncaughtAtOpen_;

    bool closed = true;
    if (bgzfO_ != nullptr) {
        closed = bgzf_close(bgzfO_) == 0;
    } else if (fpO.is_open()) {
        fpO.close();
        closed = !fpO.fail();
    } else if (outputFileName_.empty()) {
        closed = static_cast<bool>(cout.flush());
    }

    if ((unwinding || !closed) && !outputFileName_.empty()) {
        remove_regular_file(outputFileName_);
    }
    if (!unwinding && !closed) {
        throw write_error(outputFileName_);
    }
}

// Strip newlines at both ends of outTxt and write it as one block ending in a newline; empty text writes nothing.
void SAVE::save(
    string & outTxt
)
{
    const size_t last = outTxt.find_last_not_of('\n');
    if (last == string::npos) {
        outTxt.clear();
        return;
    }
    outTxt.erase(last + 1);
    outTxt.erase(0, outTxt.find_first_not_of('\n'));

    outTxt += "\n";
    write(outTxt);
    if (outputFileName_.empty() && !cout.flush()) {
        throw write_error(outputFileName_);
    }
}

// Write text exactly as given (FASTA/FASTQ records); throw if the write fails.
void SAVE::write(
    const string & text
)
{
    if (bgzfO_ != nullptr) {
        if (bgzf_write(bgzfO_, text.data(), text.size()) < 0) {
            throw write_error(outputFileName_);
        }
    } else if (outputFileName_.size() > 0) {
        if (!fpO.write(text.data(), text.size())) {
            throw write_error(outputFileName_);
        }
    } else if (!cout.write(text.data(), text.size())) {
        throw write_error(outputFileName_);
    }
}
