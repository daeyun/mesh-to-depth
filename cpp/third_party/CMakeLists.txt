#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(BUILD_SHARED_LIBS OFF)

find_package(Threads REQUIRED)
#find_package(Boost REQUIRED COMPONENTS system filesystem)


add_library(eigen3 INTERFACE)
target_include_directories(eigen3 BEFORE INTERFACE
        ${CMAKE_SOURCE_DIR}/third_party/repos/eigen
        )

add_library(cpp_gsl INTERFACE)
target_include_directories(cpp_gsl BEFORE INTERFACE
        ${CMAKE_SOURCE_DIR}/third_party/repos/GSL/include
        )

add_library(spdlog INTERFACE)
target_include_directories(spdlog BEFORE INTERFACE
        ${CMAKE_SOURCE_DIR}/third_party/repos/spdlog/include
        )
target_link_libraries(spdlog INTERFACE Threads::Threads)

add_library(nanort INTERFACE)
target_include_directories(nanort BEFORE INTERFACE
        ${CMAKE_SOURCE_DIR}/third_party/repos/nanort
        )
