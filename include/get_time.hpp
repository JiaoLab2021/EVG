#pragma once

#include <ctime>
#include <string>

// Current local time as "YYYY-mm-dd HH:MM:SS".
inline std::string getTime() {
    std::time_t now = std::time(nullptr);
    char text[64];
    std::strftime(text, sizeof(text), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    return text;
}
