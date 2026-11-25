#pragma once
#include <algorithm>
#include <vector>
#include <thread>
#include <future>
#include <exception>
#include <functional>
#include <atomic>
#include "logger.hpp"
#include <type_traits>
#include <utility>
#include <stdexcept>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include "BlockingQueue.hpp"

class ThreadPool {
public:
    explicit ThreadPool(size_t n_threads = std::thread::hardware_concurrency(), size_t queue_capacity = 0)
        : tasks_(compute_capacity_(n_threads, queue_capacity))
    {
        if (n_threads == 0) n_threads = 1;
        n_threads_ = n_threads;
        for (size_t i = 0; i < n_threads_; ++i) {
            workers_.emplace_back([this]{ this->worker_loop_(); });
        }
    }

    ~ThreadPool() {
        stop();
    }

    ThreadPool(const ThreadPool&)            = delete;
    ThreadPool(ThreadPool&&)                 = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool&&)      = delete;

    template<class F, class... Args>
    auto submit(F&& f, Args&&... args)
        -> std::future<std::invoke_result_t<F,Args...>>
    {
        using R = std::invoke_result_t<F,Args...>;
        if (stopped_.load(std::memory_order_acquire)) {
            throw std::runtime_error("submit on stopped ThreadPool");
        }

        auto task = std::make_shared<std::packaged_task<R()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        std::future<R> fut = task->get_future();

        tasks_.push([task]{ (*task)(); });

        return fut;
    }

    template<class F>
    void parallel_for(size_t count, F&& function) {
        if (count == 0) return;

        std::atomic<size_t> next{0};
        const size_t worker_count = std::min(count, n_threads_);
        std::vector<std::future<void>> futures;
        futures.reserve(worker_count);

        for (size_t worker = 0; worker < worker_count; ++worker) {
            futures.emplace_back(submit([&next, count, &function] {
                while (true) {
                    const size_t index = next.fetch_add(1, std::memory_order_relaxed);
                    if (index >= count) return;
                    std::invoke(function, index);
                }
            }));
        }

        std::exception_ptr error;
        for (auto& future : futures) {
            try {
                future.get();
            } catch (...) {
                if (!error) error = std::current_exception();
            }
        }
        if (error) std::rethrow_exception(error);
    }

    void stop() {
        bool expected = false;
        if (!stopped_.compare_exchange_strong(expected, true)) {
            return; // already stopped
        }
        tasks_.shutdown();
        for (auto& t : workers_) {
            if (t.joinable()) t.join();
        }
    }

    size_t get_pending()    const { return tasks_.size(); }
    size_t thread_count()   const { return n_threads_; }
    size_t queue_capacity() const { return tasks_.capacity(); }

private:
    static size_t compute_capacity_(size_t& n_threads, size_t user_cap) {
        if (n_threads == 0) n_threads = 1;
        if (user_cap == 0) {
            // Default capacity = 2 * n_threads
            return 2 * n_threads;
        }
        return user_cap;
    }

    void worker_loop_( ) {
        while (true) {
            auto opt = tasks_.pop();
            if (!opt) break;  // queue closed & empty
            try {
                (*opt)();
            } catch (const std::exception& e) {
                error_stream() << "Uncaught exception: " << e.what() << '\n';
            } catch (...) {
                error_stream() << "Unknown exception\n";
            }
        }
    }

    std::vector<std::thread>             workers_;
    BlockingQueue<std::function<void()>> tasks_;     // bounded queue
    std::atomic<bool>                    stopped_{false};
    size_t                               n_threads_{0};
};