project (OpenGR-AccelLib)

set(accel_ROOT "gr/accelerators")

#############################################

set(accel_relative_INCLUDE
    ${accel_ROOT}/kdtree.h
    ${accel_ROOT}/pairExtraction/bruteForceFunctor.h
    ${accel_ROOT}/pairExtraction/intersectionFunctor.h
    ${accel_ROOT}/pairExtraction/intersectionNode.h
    ${accel_ROOT}/pairExtraction/intersectionPrimitive.h
    ${accel_ROOT}/normalset.h
    ${accel_ROOT}/normalset.hpp
    ${accel_ROOT}/normalset.hpp
    ${accel_ROOT}/utils.h
    ${accel_ROOT}/gret_sdp/wrappers.h
    ${accel_ROOT}/gret_sdp/SDPAWrapper.h
    ${accel_ROOT}/gret_sdp/MOSEKWrapper.h)

if(OpenGR_USE_CHEALPIX)
    include_directories(${Chealpix_INCLUDE_DIR})
    set(APPEND accel_relative_INCLUDE ${accel_ROOT}/normalHealSet.h)
endif(OpenGR_USE_CHEALPIX)

# compute accel_INCLUDE to get absolute filenames
set(accel_INCLUDE "")
set(accel_installed_INCLUDE "")
foreach(f ${accel_relative_INCLUDE})
    list(APPEND accel_INCLUDE "${SRC_DIR}/${f}")
    list(APPEND accel_installed_INCLUDE "include/${f}")
endforeach(f)

add_library(accel INTERFACE)
target_sources(accel INTERFACE
    $<BUILD_INTERFACE:${accel_INCLUDE}>
    $<INSTALL_INTERFACE:${accel_installed_INCLUDE}>
    )
add_dependencies(accel gr::utils)
set_target_properties(accel PROPERTIES
  INTERFACE_COMPILE_FEATURES cxx_std_17
)

if(OpenGR_USE_CHEALPIX)
    add_dependencies(accel chealpix)
    target_link_libraries(accel INTERFACE chealpix)
endif(OpenGR_USE_CHEALPIX)

set(targets_export_name "${PROJECT_NAME}Targets")
install( TARGETS accel
    EXPORT "${targets_export_name}"
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include  )

add_library(gr::accel ALIAS accel)

#############################################
#
# HACK: have the files showing in the IDE, under the name 'opengr_accel_IDE'
# Source: https://stackoverflow.com/questions/27039019/how-to-have-cmake-show-headers-that-are-not-part-of-any-binary-target-in-the-ide/29214327#29214327
#
add_custom_target(opengr_accel_IDE SOURCES ${accel_INCLUDE})
