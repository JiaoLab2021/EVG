cmake_minimum_required(VERSION 3.12)

project(EVG CXX)

# 设置 C++ 标准版本 (这里使用 C++14)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 添加源文件
set(SRC_FILES_graphvcf
    src/check_vcf_sort.cpp
    src/get_time.cpp
    src/save.cpp
    src/sort_find.cpp
    src/strip_split_join.cpp
    src/delly2vcf.cpp
    src/vcf_breakpoint.cpp
    src/vcf_convert.cpp
    src/vcf_count.cpp
    src/vcf_filter.cpp
    src/vcf_merge.cpp
    src/vcf_open.cpp
    src/vcf_recall.cpp
    src/vcf_split.cpp
    src/vcf_win.cpp
)

set(SRC_FILES_fastAQ
    src/get_time.cpp
    src/save.cpp
    src/sort_find.cpp
    src/strip_split_join.cpp
    src/convert.cpp
    src/count.cpp
)

# 添加可执行文件
add_executable(graphvcf graphvcf.cpp)
add_executable(fastAQ fastAQ.cpp)

# 添加源文件
target_sources(graphvcf PRIVATE ${SRC_FILES_graphvcf})
target_sources(fastAQ PRIVATE ${SRC_FILES_fastAQ})

# 为目标添加编译选项
target_compile_options(graphvcf PRIVATE -O3)
target_compile_options(fastAQ PRIVATE -O3)

# 为目标添加编译选项（包括 -O3、额外的警告和错误和CPU优化）
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # target_compile_options(graphvcf PRIVATE -O3 -Wall -Wextra -Werror -pedantic-errors -march=native)
    # target_compile_options(fastAQ PRIVATE -O3 -Wall -Wextra -Werror -pedantic-errors -march=native)
    target_compile_options(graphvcf PRIVATE -O3 -Wall -Wextra -march=native)
    target_compile_options(fastAQ PRIVATE -O3 -Wall -Wextra -march=native)
endif()


# 添加 include 目录
target_include_directories(graphvcf PRIVATE include)
target_include_directories(fastAQ PRIVATE include)

# 添加链接库
find_package(Threads REQUIRED)
target_link_libraries(graphvcf Threads::Threads)
target_link_libraries(fastAQ Threads::Threads)

find_package(ZLIB REQUIRED)
target_link_libraries(graphvcf ZLIB::ZLIB)
target_link_libraries(fastAQ ZLIB::ZLIB)

target_link_libraries(graphvcf stdc++fs)
target_link_libraries(fastAQ stdc++fs)

# install function
install(TARGETS graphvcf DESTINATION bin)
install(TARGETS fastAQ DESTINATION bin)

# 安装头文件
install(DIRECTORY include/ DESTINATION include)

install(FILES EVG.py
        DESTINATION bin
)

# copy
install(CODE "execute_process(COMMAND cp -r src ${CMAKE_INSTALL_PREFIX}/bin/)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/EVG.py)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/EVG.py ${CMAKE_INSTALL_PREFIX}/bin/EVG)")