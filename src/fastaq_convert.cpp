#include "../include/fastaq_convert.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

#include "../include/save.hpp"

// Store the input and output file names and the thread count.
FastaConvert::FastaConvert(const std::string& fastaFilename, const std::string& outputFilename, uint32_t threads)
    : fastaFilename_(fastaFilename), outputFilename_(outputFilename), threads_(threads) {}

// Join the sequence lines of each record, skip empty lines and write the text in 1 MB chunks (.gz compressed on
// threads). Bytes equal the former single save(): no newline before the first header, one newline at the end.
void FastaConvert::convert() {
    std::ifstream fastaFile;
    fastaFile.open(fastaFilename_, std::ios::in);
    if (!fastaFile.is_open()) {
        throw std::runtime_error("'" + fastaFilename_ + "': No such file or directory.");
    }

    SAVE output(outputFilename_, threads_);
    std::string informations;
    std::string text;
    bool wroteLine = false;
    bool endsWithHeader = false;
    while (std::getline(fastaFile, informations)) {
        if (informations.empty()) {
            continue;
        }
        if (informations.find(">") != std::string::npos) {
            if (wroteLine) text += '\n';
            text += informations;
            text += '\n';
            endsWithHeader = true;
        } else {
            text += informations;
            endsWithHeader = false;
        }
        wroteLine = true;
        if (text.size() >= (1 << 20)) {
            output.write(text);
            text.clear();
        }
    }
    if (wroteLine && !endsWithHeader) text += '\n';
    output.write(text);
}
