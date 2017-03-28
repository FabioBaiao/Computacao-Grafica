cmake_minimum_required(VERSION 2.8)

# Project Name
PROJECT(cg_assignment)

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIR})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})

if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tinyxml)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tinyxml)
add_definitions(-std=c++11 -g)

add_executable(engine engine.cpp)
add_executable(engine2 engine_tinyxml2.cpp tinyxml2.cpp)
add_executable(gen gen.cpp)

target_link_libraries(engine ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} tinyxml)
target_link_libraries(engine2 ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
