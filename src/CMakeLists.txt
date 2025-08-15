file(GLOB SOURCES ./*.c ./*.h)  
add_executable(kmer-compressor ${SOURCES})  

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
    message(STATUS "ZLIB found: ${ZLIB_LIBRARIES}")
    message(STATUS "ZLIB include directories: ${ZLIB_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR, "ZLIB not found")
endif()
target_link_libraries(kmer-compressor PUBLIC ZLIB::ZLIB)
target_include_directories(kmer-compressor PUBLIC ${CMAKE_CURRENT_LIST_DIR}
    ${ZLIB_INCLUDE_DIRS})
