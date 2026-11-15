cmake_minimum_required(VERSION 3.16...3.27)
project(CosmicMuonSim CXX)

find_package(Geant4 REQUIRED)

# Generator sources live at the repository root
set(MEGEN_DIR ${PROJECT_SOURCE_DIR}/../.. CACHE PATH "MuonElectronGenerator source directory")
set(MEGEN_SOURCES ${MEGEN_DIR}/MuonElectronGenerator.cc ${MEGEN_DIR}/MT_Random.cc)
set_source_files_properties(${MEGEN_SOURCES} PROPERTIES COMPILE_OPTIONS "-w")  # legacy code, quiet build

file(GLOB APP_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
add_executable(cosmic_muons cosmic_muons.cc ${APP_SOURCES} ${MEGEN_SOURCES})
target_include_directories(cosmic_muons PRIVATE ${PROJECT_SOURCE_DIR}/include ${MEGEN_DIR})
target_link_libraries(cosmic_muons PRIVATE ${Geant4_LIBRARIES})
