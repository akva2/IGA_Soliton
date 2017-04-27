project(Kuramoto_Sivashinsky)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

find_package(TinyXML REQUIRED)
find_package(Petsc REQUIRED)
include_directories(${TINYXML_INCLUDE_DIRS} ${PETSC_INCLUDE_DIRS})

enable_language(CXX)
cmake_minimum_required(VERSION 3.5)

set(IGA_1D  IGA_1D/1_Spline_Code.cpp)

set(KURA_SIVA  KuramotoSivashinsky/Main_File.cpp)

set(EXPREVAL_SRCS expreval/except.cpp
                  expreval/expr.cpp
                  expreval/func.cpp
                  expreval/funclist.cpp
                  expreval/node.cpp
                  expreval/parser.cpp
                  expreval/vallist.cpp)
add_library(expreval ${EXPREVAL_SRCS})

add_executable(KuramotoSivashinsky ${IGA_1D} ${KURA_SIVA})
target_link_libraries(KuramotoSivashinsky expreval ${TINYXML_LIBRARIES} ${PETSC_LIBRARIES})
