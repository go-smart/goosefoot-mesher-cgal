# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()
set(External_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
message(STATUS "Pleae make sure that you have compiled Protobuf and copied libprotobuf.lib, 
protoc.exe, libprobuf-lite.lib, libprotoc.lib to 
${External_DIR}/include/arch/${CMAKE_BUILD_TYPE}/ folder, where arch=Win32 or x64")

if(MSVC)
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		if(WIN32)
			set(PROTOBUF_LIBRARY ${External_DIR}/include/Win32/Debug/libprotobuf.lib)
			set(PROTOBUF_PROTOC_EXECUTABLE ${External_DIR}/include/Win32/Debug/protoc.exe)
			set(PROTOBUF_LITE_LIBRARY ${External_DIR}/include/Win32/Debug/libprotobuf-lite.lib)
			set(PROTOBUF_PROTOC_LIBRARY ${External_DIR}/include/Win32/Debug/libprotoc.lib)
		else(WIN32)
			set(PROTOBUF_LIBRARY ${External_DIR}/include/x64/Debug/libprotobuf.lib)
			set(PROTOBUF_PROTOC_EXECUTABLE ${External_DIR}/include/x64/Debug/protoc.exe)
			set(PROTOBUF_LITE_LIBRARY ${External_DIR}/include/x64/Debug/libprotobuf-lite.lib)
			set(PROTOBUF_PROTOC_LIBRARY ${External_DIR}/include/x64/Debug/libprotoc.lib)
		endif(WIN32)
		set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
		set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS  "Debug" "Release")
	endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
	if(CMAKE_BUILD_TYPE STREQUAL "Release")
		if(WIN32)
			set(PROTOBUF_LIBRARY ${External_DIR}/include/Win32/Release/libprotobuf.lib)
			set(PROTOBUF_PROTOC_EXECUTABLE ${External_DIR}/include/Win32/Release/protoc.exe)
			set(PROTOBUF_LITE_LIBRARY ${External_DIR}/include/Win32/Release/libprotobuf-lite.lib)
			set(PROTOBUF_PROTOC_LIBRARY ${External_DIR}/include/Win32/Release/libprotoc.lib)
		else(WIN32)
			set(PROTOBUF_LIBRARY ${External_DIR}/include/x64/Release/libprotobuf.lib)
			set(PROTOBUF_PROTOC_EXECUTABLE ${External_DIR}/include/x64/Release/protoc.exe)
			set(PROTOBUF_LITE_LIBRARY ${External_DIR}/include/x64/Release/libprotobuf-lite.lib)
			set(PROTOBUF_PROTOC_LIBRARY ${External_DIR}/include/interface/x64/Release/libprotoc.lib)
		endif(WIN32)
		set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
		set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")
	endif(CMAKE_BUILD_TYPE STREQUAL "Release")

	set(PROTOBUF_INCLUDE_DIR ${External_DIR}/interface/msvc)
	set(PROTOBUF_SRC_ROOT_FOLDER ${External_DIR}/include/src)
	if(WIN32)
		include_directories(C:/Qt/4.8.6/Win32/include)
		SET(VTK_DIR C:/VTK32/bin CACHE PATH "VTK directory override" FORCE)
	else(WIN32)
		include_directories(C:/Qt/4.8.6/x64/include)
		SET(VTK_DIR C:/VTK64/bin CACHE PATH "VTK directory override" FORCE)
	endif(WIN32)
	include_directories(${External_DIR}/include/msvc ${External_DIR}/include/msvc/include)
elseif (CMAKE_COMPILER_IS_GNUCXX)
	include_directories(include/gnu)
	set(PROTOBUF_INCLUDE_DIR include/gnu)
	message("Using g++, turning on all warnings")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
endif()

find_package(Protobuf REQUIRED)
PROTOBUF_GENERATE_CPP(CGAL_SETTINGS_PROTO_SRC CGAL_SETTINGS_PROTO_HEADER cgalsettings.proto)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
set(CGAL_SETTINGS_PROTO_HEADER_DIR ${CMAKE_CURRENT_BINARY_DIR}/../MeshGenCGAL)
include_directories(${CGAL_SETTINGS_PROTO_HEADER_DIR})


#set ( Boost_USE_STATIC_LIBS OFF )
#set ( Boost_USE_MULTITHREADED ON )
#set ( Boost_USE_STATIC_RUNTIME OFF )
#include_directories(../../include)
#include_directories(../../../Triangulation_3/include)
#include_directories(../../../STL_Extension/include)
#include_directories(../../../AABB_tree/include)
add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_VERBOSE
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if ( MESH_3_VERBOSE )
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()
find_package(CGAL COMPONENTS ImageIO)
include( ${CGAL_USE_FILE} )


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})


include_directories (BEFORE include)
if (MSVC)
    find_package(Boost)
    IF (Boost_FOUND)
        INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
        ADD_DEFINITIONS( "-DHAS_BOOST" )
    ENDIF()
	include_directories(../interface/msvc/include)
	SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "C:/Boost")
	SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "C:/Boost/lib32-msvc")
else (MSVC)
    find_package(Boost COMPONENTS thread program_options filesystem REQUIRED)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
endif(MSVC)


add_library(numa-mesher-cgal STATIC
    mesher_cgal.cpp
    #visualization.cpp
    simplify.cpp
    ${CGAL_SETTINGS_PROTO_SRC}
)
target_link_libraries(numa-mesher-cgal ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${PROTOBUF_LIBRARY})

if(VTK_LIBRARIES)
  target_link_libraries(numa-mesher-cgal ${VTK_LIBRARIES})
else()
  target_link_libraries(numa-mesher-cgal vtkHybrid)
endif()

add_executable(mesher_cgal.pb main_pb.cpp)
target_link_libraries(mesher_cgal.pb numa-mesher-cgal)

if (NOT MSVC)
    add_executable(mesher_cgal main.cpp)
    target_link_libraries(mesher_cgal numa-mesher-cgal boost_program_options boost_filesystem)
    INSTALL(TARGETS mesher_cgal RUNTIME DESTINATION "bin" LIBRARY DESTINATION "lib") 
endif()


if(MSVC)
	INSTALL(TARGETS mesher_cgal.pb RUNTIME DESTINATION "bin" ARCHIVE DESTINATION "lib") 
#	INSTALL(FILES mesher_cgal.pb.exe DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#	INSTALL(FILES libmesher_cgal.lib DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()
