cmake_minimum_required(VERSION 3.15)

project(KmerCompression VERSION 0.0.1 DESCRIPTION "Test of k-mer compression and counting methods" LANGUAGES C)






set (TARGET_GROUP production CACHE STRING "What Should I Build?")  # Modules that exist in both our release and our unit tests, are added now   

if(TARGET_GROUP STREQUAL release) # If we had any files that ONLY were compiled into our release, for example startup code,   
    # it's really handy to put it here.  
    set(CMAKE_C_FLAGS "-Wall -O3")
    
    add_subdirectory(src)
 
elseif(TARGET_GROUP STREQUAL debug)    
    set(CMAKE_C_FLAGS "-Wall -Wextra -Werror -fanalyzer -g") 
    add_subdirectory(src)

elseif(TARGET_GROUP STREQUAL test) # Our tests need to let CTest know we're running tests    
    set(CMAKE_C_FLAGS "-Wall -g")
    add_definitions(-DDEBUG)
    include(CTest)  # It also needs to understand we have code in source and test directories  
    add_subdirectory(Unity)    
    add_subdirectory(src)    
    add_subdirectory(test)  
else()    
    message(FATAL_ERROR "I don't know the TARGET_GROUP you gave me!") 
endif()
