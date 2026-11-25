#pragma once

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <iomanip>
#include <unordered_map>

#include "strip_split_join.hpp"
#include "logger.hpp"

using namespace std;

struct BGZF;

// Write text to a BGZF file (name contains .gz/.GZ, compression level 5), a plain file, or stdout (empty name).
// BGZF bytes depend only on the text written, never on how it was split into calls or on the thread count.
// A SAVE destroyed while an exception propagates removes its file (regular files only), so a failed run leaves no partial output.
class SAVE
{
private:
    string outputFileName_;

    ofstream fpO;
    BGZF* bgzfO_;
    int uncaughtAtOpen_;  // exceptions already propagating when the file was opened
public:
    SAVE(
        const string & outputFileName,
        uint32_t threads = 1
    );

    ~SAVE() noexcept(false);

    SAVE(const SAVE&) = delete;
    SAVE& operator=(const SAVE&) = delete;

    void save(string & outTxt);

    void write(const string & text);
};
