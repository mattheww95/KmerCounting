file(GLOB SOURCES ./*.c)  
add_executable(kmer-compressor ${SOURCES})  
target_include_directories(kmer-compressor PUBLIC ${CMAKE_CURRENT_LIST_DIR})
