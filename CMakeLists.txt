cmake_minimum_required(VERSION 3.12)
project(EVG CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加所有的库
set(LIBRARIES check_vcf_sort strip_split_join get_time save sort_find
              vcf_breakpoint vcf_convert vcf_count vcf_filter vcf_merge
              vcf_open vcf_recall vcf_split vcf_win delly2vcf count convert)

foreach(library IN LISTS LIBRARIES)
    add_library(${library} SHARED src/${library}.cpp)
    target_compile_features(${library} PRIVATE cxx_std_14)
    target_compile_options(${library} PRIVATE -Wall -Wextra -O3)
endforeach()

# 添加可执行文件
add_executable(graphvcf graphvcf.cpp)
add_executable(fastAQ fastAQ.cpp)

# 添加头文件搜索路径
target_include_directories(graphvcf PRIVATE include)
target_include_directories(fastAQ PRIVATE include)

# 查找线程库和zlib库
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# 设置公共链接库列表
set(COMMON_LIBRARIES 
    strip_split_join
    get_time
    Threads::Threads
    ZLIB::ZLIB
    stdc++fs
)

# 添加链接库
target_link_libraries(graphvcf PUBLIC 
    check_vcf_sort
    save
    sort_find
    delly2vcf
    vcf_breakpoint
    vcf_convert
    vcf_count
    vcf_filter
    vcf_merge
    vcf_open
    vcf_recall
    vcf_split
    vcf_win
    ${COMMON_LIBRARIES}	
)

target_link_libraries(fastAQ PUBLIC 
    convert
    count
    ${COMMON_LIBRARIES}
)

# 添加编译选项
target_compile_options(graphvcf PRIVATE -O3)
target_compile_options(fastAQ PRIVATE -O3)

# 添加相对的动态库搜索路径
target_link_directories(graphvcf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_directories(fastAQ PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# 添加安装命令
install(TARGETS graphvcf fastAQ ${LIBRARIES}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
)

install(FILES EVG.py
        DESTINATION bin
)

# copy
install(CODE "execute_process(COMMAND cp -r src ${CMAKE_INSTALL_PREFIX}/bin/)")
install(CODE "execute_process(COMMAND chmod +x ${CMAKE_INSTALL_PREFIX}/bin/EVG.py)")
install(CODE "execute_process(COMMAND ln -s ${CMAKE_INSTALL_PREFIX}/bin/EVG.py ${CMAKE_INSTALL_PREFIX}/bin/EVG)")
