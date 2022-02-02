project (OpenGR-AlgoLib)

set(algo_files  AbstractMatchBase.h
                Functor4pcs.h
                FunctorBrute4pcs.h
                FunctorSuper4pcs.h
                NAryMatchBase.h
                PairwiseMatchBase.h
                PairwiseMatchBase.hpp
                PointPairFilter.h
                congruentSetExplorationBase.h
                congruentSetExplorationBase.hpp
                match3pcs.h
                match3pcs.hpp
                match4pcsBase.h
                match4pcsBase.hpp
                pairCreationFunctor.h
        )

if(OpenGR_USE_EXPERIMENTAL_FEATURES)
    set(APPEND algo_files
            experimental/GRET_SDP.h
            experimental/GRET_SDP.hpp
            )
endif()

# compute accel_INCLUDE to get absolute filenames
set(algo_INCLUDE "")
set(algo_installed_INCLUDE "")
foreach(f ${algo_files})
    list(APPEND algo_installed_INCLUDE "include/gr/algorithms/${f}")
    list(APPEND algo_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
endforeach(f)

add_library(algo INTERFACE)
target_sources(algo INTERFACE
    $<BUILD_INTERFACE:${algo_INCLUDE}>
    $<INSTALL_INTERFACE:${algo_installed_INCLUDE}>
    )
add_dependencies(algo gr::accel gr::utils)
set_target_properties(algo PROPERTIES
  INTERFACE_COMPILE_FEATURES cxx_std_17
)

set(targets_export_name "${PROJECT_NAME}Targets")
install( TARGETS algo
    EXPORT "${targets_export_name}"
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include  )

add_library(gr::algo ALIAS algo)

#############################################
#
# HACK: have the files showing in the IDE, under the name 'opengr_algo_IDE'
# Source: https://stackoverflow.com/questions/27039019/how-to-have-cmake-show-headers-that-are-not-part-of-any-binary-target-in-the-ide/29214327#29214327
#
add_custom_target(opengr_algo_IDE SOURCES ${algo_INCLUDE})
