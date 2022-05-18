# $Id: CMakeLists.txt 94031 2015-11-05 11:54:38Z ihrivnac $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B3a)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleB3a exampleB3a.cc ${sources} ${headers})
target_link_libraries(exampleB3a ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B3a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
# set(EXAMPLEB3_SCRIPTS
#   debug.mac
#   exampleB3.in
#   exampleB3.out
#   init_vis.mac
#   run1.mac
#   run2.mac
#   vis.mac
#   )

# file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/*.mac 
#                     ${PROJECT_SOURCE_DIR}/*.in
#                     ${PROJECT_SOURCE_DIR}/*.C 
# 		    ${PROJECT_SOURCE_DIR}/*.py # copying the python file to the compling folder as well
# 		    ${PROJECT_SOURCE_DIR}/*.txt # 
# 	            ${PROJECT_SOURCE_DIR}/*.mat
#                     ${PROJECT_SOURCE_DIR}/README
#                     )
# 
# foreach(_script ${EXAMPLEB3_SCRIPTS})
#   configure_file(
#     ${PROJECT_SOURCE_DIR}/${_script}
#     ${PROJECT_BINARY_DIR}/${_script}
#     COPYONLY
#     )
# endforeach()


file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/*.mac 
                    ${PROJECT_SOURCE_DIR}/*.in
                    ${PROJECT_SOURCE_DIR}/*.C 
		    ${PROJECT_SOURCE_DIR}/*.py # copying the python file to the compling folder as well
		    ${PROJECT_SOURCE_DIR}/*.txt # 
	            ${PROJECT_SOURCE_DIR}/*.mat
		    ${PROJECT_SOURCE_DIR}/*.sh
                    ${PROJECT_SOURCE_DIR}/README
                    )

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build DICOM. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

# # the macros
# set(DICOM_MACROS
#     run.mac vis.mac
#     )


# original set of DICOM data
    set(DICOM_SCRIPTS
     1.dcm 2.dcm 3.dcm 
     1.g4  2.g4 3.g4 
     1.g4dcm 2.g4dcm 3.g4dcm
     ColourMap.dat CT2Density.dat 
   )
 
# new DICOM data (in share directory)
set(DICOM_SHARE
    AltData.dat SixSlice.dat
    IM-0003-0001.dcm IM-0003-0003.dcm IM-0003-0005.dcm IM-0003-0007.dcm IM-0003-0009.dcm
    IM-0003-0002.dcm IM-0003-0004.dcm IM-0003-0006.dcm IM-0003-0008.dcm IM-0003-0010.dcm
  )

foreach(_script ${DICOM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

foreach(_script ${DICOM_SHARE})
  configure_file(
    ${PROJECT_SOURCE_DIR}/share/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

# foreach(_script ${DICOM_MACROS})
#   configure_file(
#     ${PROJECT_SOURCE_DIR}/${_script}
#     ${PROJECT_BINARY_DIR}/${_script}
#     COPYONLY
#     )
# endforeach()

if("$ENV{DICOM_USE_DCMTK}")

set_source_files_properties(
     ${sources}
     PROPERTIES COMPILE_DEFINITIONS G4_DCMTK 

)
endif()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(B3a DEPENDS exampleB3a)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleB3a DESTINATION bin )
