cmake_minimum_required(VERSION 3.16...3.27)
project(MPPCSim)

#-------------------------------------------------------------------------------
# Disable G4MT (Multi-threading support)
set(GEANT4_USE_MT OFF)
add_definitions(-DG4NOTHREADS)

#-------------------------------------------------------------------------------
# Enable GDML support
add_compile_definitions(GEANT4_USE_GDML)
option(GEANT4_USE_GDML "Use GDML Option" ON)
message(STATUS "GEANT4_USE_GDML: ${GEANT4_USE_GDML}")

#-------------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#-------------------------------------------------------------------------------
# Find ROOT package
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

#-------------------------------------------------------------------------------
# Setup include directories
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Locate sources and headers
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# If GDML is enabled, add necessary files
if(GEANT4_USE_GDML)
  file(GLOB detectors ${PROJECT_SOURCE_DIR}/gdml/*.gdml)
endif()

file(GLOB macros ${PROJECT_SOURCE_DIR}/*.mac)
file(GLOB inputs ${PROJECT_SOURCE_DIR}/*.in)

#-------------------------------------------------------------------------------
# Add the executable and link it to the necessary libraries
add_executable(MPPCSim main.cc ${sources} ${headers})
target_link_libraries(MPPCSim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#-------------------------------------------------------------------------------
# Copy necessary scripts to the build directory
set(MPPCSim_SCRIPTS
    optPhoton.mac
    gui.mac
    vis.mac
)

foreach(_script ${MPPCSim_SCRIPTS})
  configure_file(${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

#-------------------------------------------------------------------------------
# Install the executable and scripts
install(TARGETS MPPCSim DESTINATION bin)

if (GEANT4_USE_GDML)
  install(FILES ${detectors} ${macros} ${inputs} DESTINATION bin)
else()
  install(FILES ${macros} ${inputs} DESTINATION bin)
endif()
