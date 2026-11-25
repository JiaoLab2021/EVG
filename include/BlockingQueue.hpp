#pragma once
/*-----------------------------------------------------------------------
 * BlockingQueue<T>
 * A thread-safe FIFO queue with shutdown support and optional capacity.
 *
 * - If capacity_ == 0, the queue is unbounded (push never blocks).
 * - If capacity_ > 0, push() blocks when the queue is full until an item
 *   is popped or shutdown() is called.
 *
 * 2025-10-24: Added capacity + "not full" condition variable to support
 * bounded producer behavior (back-pressure).
 *---------------------------------------------------------------------*/
#include <queue>
#include <mutex>
#include <condition_variable>
#include <optional>
#include <stdexcept>
#include <cstddef>

template<class T>
class BlockingQueue {
public:
    // Construct with an optional capacity; 0 means unbounded.
    explicit BlockingQueue(size_t capacity = 0)
        : capacity_(capacity) {}

    BlockingQueue(const BlockingQueue&)            = delete;
    BlockingQueue& operator=(const BlockingQueue&) = delete;

    /*-------------------------------------------------------------------
     * push(U&& value)
     * - Blocks if the queue is full (when capacity_ > 0).
     * - Throws std::runtime_error if shutdown() has been called.
     *------------------------------------------------------------------*/
    template<class U>
    void push(U&& value) {
        std::unique_lock<std::mutex> lock(mtx_);
        // Wait until: (shutdown) OR (unbounded) OR (size < capacity)
        cv_not_full_.wait(lock, [this]{
            return done_ || capacity_ == 0 || q_.size() < capacity_;
        });
        if (done_) {
            throw std::runtime_error("push on stopped BlockingQueue");
        }
        q_.emplace(std::forward<U>(value));
        lock.unlock();
        cv_not_empty_.notify_one();
    }

    /*-------------------------------------------------------------------
     * try_push(U&& value)
     * - Non-blocking; returns false if full or shutdown.
     *------------------------------------------------------------------*/
    template<class U>
    bool try_push(U&& value) {
        std::lock_guard<std::mutex> lock(mtx_);
        if (done_) return false;
        if (capacity_ != 0 && q_.size() >= capacity_) return false;
        q_.emplace(std::forward<U>(value));
        cv_not_empty_.notify_one();
        return true;
    }

    /*-------------------------------------------------------------------
     * pop()
     * - Blocks until an item is available or the queue is shutdown AND
     *   empty; in the latter case returns std::nullopt.
     *------------------------------------------------------------------*/
    std::optional<T> pop() {
        std::unique_lock<std::mutex> lock(mtx_);
        cv_not_empty_.wait(lock, [this]{
            return done_ || !q_.empty();
        });
        if (q_.empty()) {
            // done_ && empty => graceful end
            return std::nullopt;
        }
        T val = std::move(q_.front());
        q_.pop();
        lock.unlock();
        cv_not_full_.notify_one();  // we consumed one slot; wake up producers
        return val;
    }

    /*-------------------------------------------------------------------
     * try_pop(T& out)
     * - Non-blocking; returns false if the queue is empty.
     *------------------------------------------------------------------*/
    bool try_pop(T& out) {
        std::lock_guard<std::mutex> lock(mtx_);
        if (q_.empty()) return false;
        out = std::move(q_.front());
        q_.pop();
        cv_not_full_.notify_one();
        return true;
    }

    /*-------------------------------------------------------------------
     * shutdown()
     * - Marks the queue as done; wakes all waiting producers/consumers.
     *------------------------------------------------------------------*/
    void shutdown() {
        {
            std::lock_guard<std::mutex> lock(mtx_);
            done_ = true;
        }
        cv_not_empty_.notify_all();
        cv_not_full_.notify_all();
    }

    /*-------------------------------------------------------------------
     * size(), empty(), capacity()
     * - Size/capacity queries (thread-safe).
     *------------------------------------------------------------------*/
    size_t size() const {
        std::lock_guard<std::mutex> lock(mtx_);
        return q_.size();
    }

    bool empty() const {
        return size() == 0;
    }

    size_t capacity() const {
        return capacity_;
    }

private:
    mutable std::mutex      mtx_;
    std::condition_variable cv_not_empty_;
    std::condition_variable cv_not_full_;
    std::queue<T>           q_;
    size_t                  capacity_ = 0;  // 0 = unbounded
    bool                    done_     = false;
};