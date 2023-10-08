# ##############################################################################
# g++ compiler configuration
# ##############################################################################

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++17")

# Add the basic compiler options for debug version
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3")
# Add the basic compiler options for release version
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -Ofast -fstrict-aliasing -march=native -DNDEBUG")