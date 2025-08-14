
add_executable(test_sequence_parsing test_sequence.c)   



find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
    message(STATUS "ZLIB found: ${ZLIB_LIBRARIES}")
    message(STATUS "ZLIB include directories: ${ZLIB_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR, "ZLIB not found")
endif()

file(GLOB SOURCES ./*.txt ./*.fasta ./*.fastq)

target_link_libraries(test_sequence_parsing Unity ZLIB::ZLIB)   
target_include_directories(test_sequence_parsing PUBLIC
    ${PROJECT_SOURCE_DIR}/src ${ZLIB_INCLUDE_DIRS})
file(COPY ${SOURCES} test.gzip DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(sequence_parsing_test test_sequence_parsing)
