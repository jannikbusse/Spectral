cmake_minimum_required(VERSION "3.16.3")

set(CMAKE_CXX_STANDARD 17)
project("spectral")
include_directories(src)
include_directories(include)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

add_executable("${PROJECT_NAME}"
    "src/main.cpp"
    )


add_definitions(-DPROJ_DIR=\"${PROJECT_SOURCE_DIR}\")
#install(TARGETS "${PROJECT_NAME}" DESTINATION bin)
target_link_libraries("${PROJECT_NAME}"  Threads::Threads)
target_link_libraries("${PROJECT_NAME}" OpenMP::OpenMP_CXX)
target_link_libraries("${PROJECT_NAME}" fftw3)
