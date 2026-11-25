#pragma once

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

class HelpPrinter {
public:
    void section(const std::string& title) {
        rows_.push_back({title + ":", "", "", true});
    }

    void line(const std::string& option, const std::string& type, const std::string& description) {
        rows_.push_back({option, type, description, false});
    }

    void blank() {
        rows_.push_back({"", "", "", true});
    }

    void print(std::ostream& out) const {
        std::size_t option_width = 0;
        std::size_t type_width = 0;
        for (const auto& row : rows_) {
            if (row.plain) continue;
            option_width = std::max(option_width, row.option.size());
            type_width = std::max(type_width, row.type.size());
        }
        option_width += 4;
        if (type_width != 0) type_width += 4;

        for (const auto& row : rows_) {
            if (row.plain) {
                out << row.option << '\n';
            } else {
                out << "  " << row.option << std::string(option_width - row.option.size(), ' ')
                    << row.type << std::string(type_width - row.type.size(), ' ')
                    << row.description << '\n';
            }
        }
    }

private:
    struct Row {
        std::string option;
        std::string type;
        std::string description;
        bool plain;
    };

    std::vector<Row> rows_;
};
