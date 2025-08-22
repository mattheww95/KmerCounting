file(GLOB SOURCES ./*.c)  

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
    message(STATUS "ZLIB found: ${ZLIB_LIBRARIES}")
    message(STATUS "ZLIB include directories: ${ZLIB_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR, "ZLIB not found")
endif()

if(TARGET_GROUP STREQUAL test)

    add_definitions(-DDEBUG)
    add_library(Testing STATIC ${SOURCES})
    target_include_directories(Testing PUBLIC ${CMAKE_CURRENT_LIST_DIR}
        ${ZLIB_INCLUDE_DIRS})
    target_link_libraries(Testing PUBLIC ZLIB::ZLIB)

elseif(TARGET_GROUP STREQUAL debug)

    add_executable(kmer-compressor-debug ${SOURCES})  
    target_include_directories(kmer-compressor-debug PUBLIC ${CMAKE_CURRENT_LIST_DIR}
        ${ZLIB_INCLUDE_DIRS})
    target_link_libraries(kmer-compressor-debug PUBLIC ZLIB::ZLIB)

elseif(TARGET_GROUP STREQUAL release)
    add_executable(kmer-compressor ${SOURCES})  
    target_link_libraries(kmer-compressor PUBLIC ZLIB::ZLIB)
    target_include_directories(kmer-compressor PUBLIC ${CMAKE_CURRENT_LIST_DIR}
        ${ZLIB_INCLUDE_DIRS})
endif()
