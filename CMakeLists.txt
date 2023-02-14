cmake_minimum_required(VERSION 3.23)
project(lognormal)

set(CMAKE_CXX_STANDARD 20)

add_executable(lognormal main.cpp)

find_package(Boost REQUIRED COMPONENTS iostreams)
target_link_libraries(lognormal Boost::iostreams)