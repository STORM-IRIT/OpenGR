cmake_minimum_required (VERSION 3.3)
project (OpenGR-Apps LANGUAGES CXX)

add_subdirectory(io)


if(OpenGR_COMPILE_APPS)
    add_subdirectory(Super4PCS)
    add_subdirectory(PCLWrapper)
    add_subdirectory(ExtPointBinding)
endif(OpenGR_COMPILE_APPS)
