#pragma once

#include <cstdint>
#include <string>

// Write a random fraction of the records of a FASTA/FASTQ file, or `frac` records when frac >= 1, to outputFileName
// (BGZF when the name contains .gz) or to stdout when it is empty; throw on any error.
void stk_sample(const std::string& inputFileName, const std::string& outputFileName, double frac, const long int& seedNum,
                const std::string& mode, uint32_t threads);
