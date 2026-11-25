#pragma once

#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "logger.hpp"
#include "task_pool.hpp"
#include "vcf_open.hpp"

// The one place of the VCF order rule. Used alone by check_vcf_sort, or inside an index pass so the file is read once.
// Everything this thread prints while the pass lives is collected into the pass's Outcome: after a sorted file it keeps all of it
// and the index's first error; for an unsorted or unreadable file only what was printed before the file was read, with that error.
// Messages, errors and exit codes therefore match "check_vcf_sort, then index".
class SortedIndexPass {
public:
    // Start collecting this thread's diagnostics and log the order check.
    explicit SortedIndexPass(const std::string& file) : file_(file), capture_(std::make_unique<DiagnosticCapture>()) {
        log_stream() << "Checking sort: '" << file_ << "'\n";
        mark_ = capture_->size();
    }

    // Throw the unsorted-file error if this data record (no '#') comes before the previous one on its chromosome.
    void check_order(const VCFINFOSTRUCT& record) {
        if (record.CHROM != preChromosome_) {
            preChromosome_ = record.CHROM;
            preRefStart_ = 0;
        }
        if (record.POS >= preRefStart_) {
            preRefStart_ = record.POS;
        } else {
            throw std::runtime_error("unsorted file -> " + file_ + " " + record.CHROM + " "
                + std::to_string(preRefStart_) + ">" + std::to_string(record.POS));
        }
    }

    // Run the index step for one record unless an earlier step failed; keep the first exception.
    template <class F>
    void index(F&& step) {
        if (error_) return;
        try {
            step();
        } catch (...) {
            error_ = std::current_exception();
        }
    }

    // The sorted file is done: value with everything printed and the first index error.
    template <class T>
    Outcome<T> finish(T value) {
        Outcome<T> outcome;
        outcome.value = std::move(value);
        outcome.diagnostics = capture_->take();
        outcome.error = error_;
        capture_.reset();
        return outcome;
    }

    // The file is unsorted or could not be read: error with only what was printed before the file was read.
    template <class T>
    Outcome<T> fail(std::exception_ptr error) {
        Outcome<T> outcome;
        outcome.diagnostics = capture_->take();
        outcome.diagnostics.resize(mark_);
        outcome.error = std::move(error);
        capture_.reset();
        return outcome;
    }

private:
    std::string file_;
    std::unique_ptr<DiagnosticCapture> capture_;
    size_t mark_ = 0;
    std::string preChromosome_;
    uint32_t preRefStart_ = 0;
    std::exception_ptr error_;
};
