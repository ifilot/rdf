# RDF (c) by Ivo Filot <i.a.w.filot@tue.nl>
#
# RDF is licensed under a
# Creative Commons Attribution 4.0 International License.
#
# You should have received a copy of the license along with this
# work. If not, see <http://creativecommons.org/licenses/by/4.0/>.
#
# set minimum cmake requirements
cmake_minimum_required(VERSION 2.8)
project (rdf)

# add custom directory to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/modules )

# set Boost
SET(BOOST_INCLUDEDIR "/usr/include")
SET(BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu")
set (Boost_NO_SYSTEM_PATHS ON)
set (Boost_USE_MULTITHREADED ON)
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_STATIC_RUNTIME OFF)
set (BOOST_ALL_DYN_LINK OFF)

# Include libraries
find_package(Boost COMPONENTS chrono regex REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(TCLAP tclap REQUIRED)

# Set include folders
include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_BINARY_DIR}
                    ${Boost_INCLUDE_DIRS}
                    ${TCLAP_INCLUDE_DIRS})

# Enable OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    option(HAS_OPENMP "OpenMP enabled" ON)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add sources
file(GLOB_RECURSE SOURCES "*.cpp")
add_executable(rdf ${SOURCES})

# Set C++14
add_definitions(-std=c++14)

# Link libraries
target_link_libraries(rdf ${Boost_LIBRARIES})

# add Boost definition
add_definitions(-DBOOST_LOG_DYN_LINK)

# add Wno-literal-suffix to suppress warning messages
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")

###
# Installing
##
install (TARGETS rdf DESTINATION bin)
