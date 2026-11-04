cmake_minimum_required(VERSION 3.13)

project(BOUT_Connect VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX20_STANDARD_COMPILE_OPTION "-std=c++20")

set(ELMPB_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/elm_pb/elm_pb.cxx)

find_package(bout++ REQUIRED)
find_package(pcms REQUIRED)


add_library(bout_connect INTERFACE)
add_library(BOUT_Connect::bout_connect ALIAS bout_connect)
target_include_directories(bout_connect INTERFACE include)
target_compile_features(bout_connect INTERFACE cxx_std_17)
target_compile_definitions(bout_connect INTERFACE BC_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")


add_library(bc_elmpb STATIC src/elmpb.cpp)
target_link_libraries(bc_elmpb PUBLIC bout_connect bout++::bout++)
target_compile_definitions(bc_elmpb PRIVATE ELMPB_SOURCE="${ELMPB_SOURCE}")


add_executable(coexist tests/coexist.cpp)
target_link_libraries(coexist PRIVATE bc_elmpb pcms::core)

add_executable(transfer tests/transfer.cpp)
target_link_libraries(transfer PRIVATE bc_elmpb pcms::core)

#add_executable(transfer_oh2d tests/transfer_oh2d.cpp)
#target_link_libraries(transfer_oh2d PRIVATE bout_connect pcms::core)
