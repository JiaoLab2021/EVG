#pragma once

#include <cstdint>
#include <string>

// Read and base totals.
struct countStruct {
    uint64_t readNum;
    uint64_t readBase;

    countStruct() : readNum(0), readBase(0) {}
};

// Count reads and bases of one FASTA/FASTQ file or a read pair and save readBase/readNum/readLen.
class Count {
private:
    std::string inputFileName1_;
    std::string inputFileName2_;
    std::string outputFileName_;
    uint32_t threads_;

    countStruct countOut_;
public:
    Count(const std::string& inputFileName1, const std::string& inputFileName2, const std::string& outputFileName,
          uint32_t threads);

    void count_reads();
    void save_result();
};
