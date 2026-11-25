#pragma once

#include <cstdint>
#include <string>

// Rewrite a FASTA file with each sequence on one line.
class FastaConvert {
private:
    std::string fastaFilename_;
    std::string outputFilename_;
    uint32_t threads_;
public:
    FastaConvert(const std::string& fastaFilename, const std::string& outputFilename, uint32_t threads);

    void convert();
};
