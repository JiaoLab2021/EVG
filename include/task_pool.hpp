#pragma once

#include <cstddef>
#include <cstdint>
#include <exception>
#include <future>
#include <string>
#include <type_traits>
#include <utility>

#include "ThreadPool.hpp"
#include "logger.hpp"

// ------ Process-wide compute pool sized from -t ------
// With the htslib pool (attach_bgzf_pool in hts_file.hpp) a process runs about 2 * threads + 2 threads: the compute
// workers, the htslib workers, the main thread and one htslib I/O thread per open threaded BGZF file.
// Pool tasks must not submit to the pool or wait for its futures (no nesting).

// Compute workers shared by every threaded step; created with `threads` workers on first use and kept until exit.
inline ThreadPool& shared_pool(uint32_t threads) {
    static ThreadPool* pool = new ThreadPool(threads);
    return *pool;
}

// Run task on the shared pool when threads > 1, otherwise right here; its result or exception comes through the future.
template <class F>
auto start_task(uint32_t threads, F&& task) -> std::future<std::invoke_result_t<F>> {
    if (threads > 1) return shared_pool(threads).submit(std::forward<F>(task));
    std::packaged_task<std::invoke_result_t<F>()> inline_task(std::forward<F>(task));
    auto future = inline_task.get_future();
    inline_task();
    return future;
}

// The result of work done on a worker (one record, one file): its value or the exception it threw, and what it printed.
template <class T>
struct Outcome {
    T value{};
    std::exception_ptr error;
    std::string diagnostics;

    // Print the diagnostics where the serial loop printed them, then rethrow the error or return the value.
    T& get() {
        diagnostic_stream() << diagnostics;
        if (error) std::rethrow_exception(error);
        return value;
    }
};
