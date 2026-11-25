#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>

// Where this thread's diagnostics go: stderr, or the buffer of an active DiagnosticCapture.
inline std::ostream*& diagnostic_target() {
    thread_local std::ostream* target = &std::cerr;
    return target;
}

inline std::ostream& diagnostic_stream() { return *diagnostic_target(); }
inline std::ostream& output_stream() { return std::cout; }

// A stream buffer that appends to a string, so the text can be measured and taken without copying.
class StringAppendBuffer : public std::streambuf {
public:
    std::string text;

protected:
    int_type overflow(int_type ch) override {
        if (!traits_type::eq_int_type(ch, traits_type::eof())) text.push_back(traits_type::to_char_type(ch));
        return traits_type::not_eof(ch);
    }

    std::streamsize xsputn(const char* s, std::streamsize n) override {
        text.append(s, static_cast<std::size_t>(n));
        return n;
    }
};

// Collect this thread's diagnostics in a buffer while the object lives, so a worker's messages can be printed later in input order.
class DiagnosticCapture {
public:
    DiagnosticCapture() : stream_(&buffer_), previous_(diagnostic_target()) { diagnostic_target() = &stream_; }
    ~DiagnosticCapture() { diagnostic_target() = previous_; }

    DiagnosticCapture(const DiagnosticCapture&) = delete;
    DiagnosticCapture& operator=(const DiagnosticCapture&) = delete;

    // The diagnostics written so far.
    const std::string& text() const { return buffer_.text; }

    // The number of characters written so far.
    std::size_t size() const { return buffer_.text.size(); }

    // Move the diagnostics written so far out of the buffer, leaving it empty.
    std::string take() {
        std::string text;
        text.swap(buffer_.text);
        return text;
    }

private:
    StringAppendBuffer buffer_;
    std::ostream stream_;
    std::ostream* previous_;
};

namespace logger_detail {

inline const auto start_wall = std::chrono::steady_clock::now();
inline const std::clock_t start_cpu = std::clock();

// Pad or cut a function name to 10 characters.
inline std::string format_function(const char* function) {
    std::string text(function);
    text.resize(10, '.');
    return text;
}

// Format a number with at least `width` characters, padding with dots.
inline std::string format_number(double number, std::size_t width) {
    std::string text = std::to_string(number);
    const auto dot = text.find('.');
    width = std::max(width, dot == std::string::npos ? text.size() : dot);
    text.resize(width, '.');
    return text;
}

// Write the "[L::function::wall*cpu] " prefix to the diagnostic stream.
inline std::ostream& write_prefix(char level, const char* function) {
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_wall).count();
    const double cpu = static_cast<double>(std::clock() - start_cpu) / CLOCKS_PER_SEC;
    const std::string prefix = std::string("[") + level + "::" + format_function(function) + "::"
        + format_number(wall, 6) + "*" + format_number(wall > 0.0 ? cpu / wall : 0.0, 4) + "] ";
    return diagnostic_stream().write(prefix.data(), static_cast<std::streamsize>(prefix.size()));
}

}

#define log_stream() logger_detail::write_prefix('M', __func__)
#define stage_stream() (diagnostic_stream().put('\n'), logger_detail::write_prefix('M', __func__))
#define warning_stream() logger_detail::write_prefix('W', __func__)
#define error_stream() logger_detail::write_prefix('E', __func__)
