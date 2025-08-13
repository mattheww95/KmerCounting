
add_executable(test_sequence_parsing test_sequence.c)   

target_link_libraries(test_sequence_parsing Unity)   

target_include_directories(test_sequence_parsing PUBLIC ${PROJECT_SOURCE_DIR}/src)

file(COPY test.gzip DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(sequence_parsing_test test_sequence_parsing)
