project (OpenGR-Tests)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${SRC_DIR})


find_package(OpenMP)
set(OpenGRTestsDeps)
if(OpenMP_CXX_FOUND)
    set(OpenGRTestsDeps OpenMP::OpenMP_CXX)
    add_definitions(-DOpenGR_USE_OPENMP)
    message(STATUS "[Tests] Parallelization enabled using OpenMP")
endif()

#############################################
## start by downloading standford datasets, that will be used
## to test the approach
include(ExternalProject)
ExternalProject_Add(dataset-standford-armadillo
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/datasets/armadillo
    URL http://graphics.stanford.edu/pub/3Dscanrep/armadillo/Armadillo_scans.tar.gz
    URL_MD5 437eafa8d3d07067e17a42c727a0c28d
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
ExternalProject_Add(dataset-standford-bunny
    SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/datasets/bunny
    URL http://graphics.stanford.edu/pub/3Dscanrep/bunny.tar.gz
    URL_MD5 232bb90cac5f2ac77140e09a59bcffe4
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

#############################################
## testing files
set (testing_SRCS
    testing.h
)

#############################################
## externalAppTestD

include(ExternalProject)
ExternalProject_Add(
        externalAppTestTarget
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/externalAppTest"

        CMAKE_ARGS
            -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
            -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
            -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
            -DCMAKE_INSTALL_MESSAGE=${EXTERNAL_INSTALL_MESSAGE}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_CXX_FLAGS=${EXTERNAL_CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}
            -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
            -DCMAKE_CXX_FLAGS_RELWITHDEBINFO=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
            -DCMAKE_SHARED_LINKER_FLAGS=${CMAKE_SHARED_LINKER_FLAGS}
            -DCMAKE_OBJECT_PATH_MAX=${CMAKE_OBJECT_PATH_MAX}
            -DCMAKE_MACOSX_RPATH=TRUE
            -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}
            -DCMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX}/lib/cmake
            -DCMAKE_INSTALL_RPATH=${CMAKE_INSTALL_PREFIX}/lib
            -DEIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}
            -DOpenGR_DIR=${CMAKE_INSTALL_PREFIX}/lib/cmake/opengr
            --no-warn-unused-cli
        # override default behaviours
        UPDATE_COMMAND ""
        EXCLUDE_FROM_ALL
)

add_dependencies(externalAppTestTarget ForceOpenGRInstallation)
add_dependencies(buildtests externalAppTestTarget)
add_test(NAME externalAppTest
         COMMAND OpenGR-externalAppTest)

#############################################
## gret-sdp
set(gret-sdp_SRCS gret-sdp.cc)
find_package(SDPA REQUIRED)
include_directories(${SDPA_INCLUDE_DIRS})

find_package(Boost 1.54 COMPONENTS system filesystem  REQUIRED)
if (Boost_FOUND)
include_directories(${Boost_INCLUDE_DIRS})

    add_executable(gret-sdp ${gret-sdp_SRCS} ${testing_SRCS})
    add_dependencies(buildtests gret-sdp)
    add_test(NAME gret-sdp
        #CONFIGURATIONS Release
        COMMAND gret-sdp)
    target_link_libraries(gret-sdp ${SDPA_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} gr::accel gr::algo gr::utils)      

endif (Boost_FOUND)

#############################################
## pair extraction
set(pair_extraction_SRCS
    pair_extraction.cc
)
add_executable(pair_extraction ${pair_extraction_SRCS} ${testing_SRCS})
add_dependencies(buildtests pair_extraction)
add_test(NAME pair_extraction
         #CONFIGURATIONS Release
         COMMAND pair_extraction)
target_link_libraries(pair_extraction gr::accel gr::algo gr::utils ${OpenGRTestsDeps})
if(OpenGR_USE_CHEALPIX)
    target_link_libraries(pair_extraction ${Chealpix_LIBS} )
endif(OpenGR_USE_CHEALPIX)

if(MSVC)
     target_compile_options(pair_extraction PRIVATE /bigobj)
 endif(MSVC)

#############################################
## quad extraction
#set(quad_extraction_SRCS
#    quad_extraction.cc
#)
#add_executable(quad_extraction ${quad_extraction_SRCS} ${testing_SRCS})
#add_dependencies(buildtests quad_extraction)
#add_test(NAME quad_extraction
#         #CONFIGURATIONS Release
#         COMMAND quad_extraction)#target_link_libraries(quad_extraction open_accel open_algo)
#if(OpenGR_USE_CHEALPIX)
#    target_link_libraries(quad_extraction ${Chealpix_LIBS} )
#endif(OpenGR_USE_CHEALPIX)


find_package(Filesystem QUIET)

#############################################
## matching
set(matching_SRCS
    matching.cc)

if (CXX_FILESYSTEM_HAVE_FS)
    add_executable(matching ${matching_SRCS} ${testing_SRCS})
    add_dependencies(buildtests matching)
    set (matching_subtests 0 1)
    foreach(subtest_id ${matching_subtests})
        add_test(NAME "matching_${subtest_id}"
                 #CONFIGURATIONS Release
                 COMMAND matching ${subtest_id})
        set_tests_properties ( "matching_${subtest_id}" PROPERTIES TIMEOUT 600) # in seconds
    endforeach()
    target_link_libraries(matching gr::accel gr::algo gr::utils gr::io ${OpenGRTestsDeps} std::filesystem)
    target_compile_definitions(matching PUBLIC -DCXX_FILESYSTEM_HAVE_FS -DCXX_FILESYSTEM_IS_EXPERIMENTAL=$<BOOL:${CXX_FILESYSTEM_IS_EXPERIMENTAL}> -DCXX_FILESYSTEM_NAMESPACE=${CXX_FILESYSTEM_NAMESPACE})
    add_dependencies(matching dataset-standford-armadillo dataset-standford-bunny)
    if(OpenGR_USE_CHEALPIX)
        target_link_libraries(matching ${Chealpix_LIBS} )
    endif(OpenGR_USE_CHEALPIX)
    if(MSVC)
        target_compile_options(matching PRIVATE /bigobj)
    endif(MSVC)
endif (CXX_FILESYSTEM_HAVE_FS)

#############################################
## matching3pcs
set(matching3pcs_SRCS
        matching3pcs.cc
        )

if (CXX_FILESYSTEM_HAVE_FS)
    add_executable(matching3pcs ${matching3pcs_SRCS} ${testing_SRCS})
    add_dependencies(buildtests matching3pcs)
    set (matching3pcs_subtests 0 1)
    foreach(subtest_id ${matching3pcs_subtests})
        add_test(NAME "matching3pcs_${subtest_id}"
                #CONFIGURATIONS Release
                COMMAND matching3pcs ${subtest_id})
        set_tests_properties ( "matching3pcs_${subtest_id}" PROPERTIES TIMEOUT 600) # in seconds
    endforeach()
    target_link_libraries(matching3pcs gr::accel gr::algo gr::utils gr::io ${OpenGRTestsDeps} std::filesystem)
    target_compile_definitions(matching3pcs PUBLIC -DCXX_FILESYSTEM_HAVE_FS -DCXX_FILESYSTEM_IS_EXPERIMENTAL=$<BOOL:${CXX_FILESYSTEM_IS_EXPERIMENTAL}> -DCXX_FILESYSTEM_NAMESPACE=${CXX_FILESYSTEM_NAMESPACE})
    add_dependencies(matching3pcs dataset-standford-armadillo dataset-standford-bunny)
    if(OpenGR_USE_CHEALPIX)
        target_link_libraries(matching3pcs ${Chealpix_LIBS} )
    endif(OpenGR_USE_CHEALPIX)
    if(MSVC)
        target_compile_options(matching3pcs PRIVATE /bigobj)
    endif(MSVC)


    #############################################
    ## io
    if (OpenGR_COMPILE_APPS)
        set(io_SRCS io.cc )
        add_executable(test_io ${io_SRCS} ${testing_SRCS})
        add_test(NAME test_io
                #CONFIGURATIONS Release
                COMMAND test_io)
        set_tests_properties ( test_io PROPERTIES TIMEOUT 600) # in seconds
        target_link_libraries(test_io gr::io ${OpenGRTestsDeps} std::filesystem)
        target_compile_definitions(test_io PUBLIC -DCXX_FILESYSTEM_HAVE_FS -DCXX_FILESYSTEM_IS_EXPERIMENTAL=$<BOOL:${CXX_FILESYSTEM_IS_EXPERIMENTAL}> -DCXX_FILESYSTEM_NAMESPACE=${CXX_FILESYSTEM_NAMESPACE})
        add_dependencies(test_io dataset-standford-bunny)
        add_dependencies(buildtests test_io)
        if(MSVC)
            target_compile_options(test_io PRIVATE /bigobj)
        endif(MSVC)
    endif (OpenGR_COMPILE_APPS)
endif (CXX_FILESYSTEM_HAVE_FS)

