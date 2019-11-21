project (OpenGR-Apps-IOLib)

# Internal static library used for apps and tests

set(io_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

file(GLOB_RECURSE io_sources ${io_ROOT}/*.cc)
file(GLOB_RECURSE io_headers ${io_ROOT}/*.h ${io_ROOT}/*.hpp)

add_library(opengr_apps_io STATIC ${io_sources} ${io_headers})
target_include_directories(opengr_apps_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )

if(IO_USE_BOOST)
    find_package(Boost 1.54 COMPONENTS system filesystem REQUIRED)
    if (Boost_FOUND)
        message(STATUS "[Deps] Use Boost ${Boost_VERSION} for image loading")
        add_definitions(-DUSE_BOOST)
        target_include_directories(opengr_apps_io PUBLIC ${Boost_INCLUDE_DIRS} )
        target_link_libraries(opengr_apps_io ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
    else ()
        message(STATUS "[Deps] Boost not found")
    endif (Boost_FOUND)
endif(IO_USE_BOOST)

add_library(gr::io ALIAS opengr_apps_io)
