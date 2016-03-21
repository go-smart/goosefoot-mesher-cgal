add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_VERBOSE
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if ( MESH_3_VERBOSE )
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()
find_package(CGAL COMPONENTS ImageIO)
include( ${CGAL_USE_FILE} )

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

find_package(Protobuf REQUIRED)
PROTOBUF_GENERATE_CPP(CGAL_SETTINGS_PROTO_SRC CGAL_SETTINGS_PROTO_HEADER ${PROTOBUF_DIR}/cgalsettings.proto)
get_filename_component(CGAL_SETTINGS_PROTO_INCLUDE_DIR ${CGAL_SETTINGS_PROTO_HEADER} DIRECTORY)

add_subdirectory(vtkIsoVolume)

message(${CGAL_SETTINGS_PROTO_SRC})
add_library(numa-mesher-cgal STATIC
    mesher_cgal.cpp
    ${CGAL_SETTINGS_PROTO_SRC}
)
target_include_directories(numa-mesher-cgal PUBLIC ${CGAL_SETTINGS_PROTO_INCLUDE_DIR})
target_link_libraries(numa-mesher-cgal ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${PROTOBUF_LIBRARY} vtkIsoVolume)

if(VTK_LIBRARIES)
  target_link_libraries(numa-mesher-cgal ${VTK_LIBRARIES})
else()
  target_link_libraries(numa-mesher-cgal vtkHybrid)
endif()
