# Source: https://github.com/coin-or/Gravity/blob/master/cmake/FindMOSEK.cmake


set(MOSEK_ROOT_DIR "$ENV{MOSEK_ROOT_DIR}" CACHE PATH "MOSEK root directory.")
message("Looking for MOSEK in ${MOSEK_ROOT_DIR}")


if(UNIX)
find_path(MOSEK_INCLUDE_DIR	fusion.h	HINTS "${MOSEK_ROOT_DIR}/9.2/tools/platform/linux64x86/h")
find_library (MOSEK_LIBRARY1 	libfusion64.so HINTS "${MOSEK_ROOT_DIR}/9.2/tools/platform/linux64x86/bin")
find_library(MOSEK_LIBRARY2 	libmosek64.so HINTS "${MOSEK_ROOT_DIR}/9.2/tools/platform/linux64x86/bin")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MOSEK DEFAULT_MSG MOSEK_LIBRARY1 MOSEK_LIBRARY2 MOSEK_INCLUDE_DIR)

if(MOSEK_FOUND)
    message("—- Found MOSEK under ${MOSEK_INCLUDE_DIR}")
    set(MOSEK_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})
    set(MOSEK_LIBRARIES ${MOSEK_LIBRARY1} ${MOSEK_LIBRARY2})
    message("—- Set MOSEK lib ${MOSEK_LIBRARY2}")
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(MOSEK_LIBRARIES "${MOSEK_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(MOSEK_FOUND)

mark_as_advanced(MOSEK_LIBRARY1 MOSEK_LIBRARY2 MOSEK_INCLUDE_DIR)
