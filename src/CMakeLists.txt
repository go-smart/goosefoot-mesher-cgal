add_executable(mesher_cgal.pb main_pb.cpp)
target_link_libraries(mesher_cgal.pb numa-mesher-cgal)

find_package(Protobuf REQUIRED)
PROTOBUF_GENERATE_CPP(CGAL_SETTINGS_PROTO_SRC CGAL_SETTINGS_PROTO_HEADER ${PROTOBUF_DIR}/cgalsettings.proto)

if (NOT MSVC)
    add_executable(mesher_cgal main.cpp ${CGAL_SETTINGS_PROTO_HEADER})
    target_include_directories(mesher_cgal PUBLIC ${CGAL_SETTINGS_PROTO_INCLUDE_DIR})
    target_link_libraries(mesher_cgal numa-mesher-cgal boost_program_options boost_filesystem)
    INSTALL(TARGETS mesher_cgal RUNTIME DESTINATION "bin" LIBRARY DESTINATION "lib") 
endif()


if(MSVC)
	INSTALL(TARGETS mesher_cgal.pb RUNTIME DESTINATION "bin" ARCHIVE DESTINATION "lib") 
#	INSTALL(FILES mesher_cgal.pb.exe DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#	INSTALL(FILES libmesher_cgal.lib DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()
