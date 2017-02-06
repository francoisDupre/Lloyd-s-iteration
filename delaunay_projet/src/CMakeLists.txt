project (delaunay)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

cmake_minimum_required(VERSION 2.4.5)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )

if( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt5)
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL)
find_package(GLUT)
find_package(OpenGL)


set( SRCS glviewer.cpp main.cpp window.cpp dt.h scene.h)


if(WIN32)
	add_definitions( "/W3 /D_CRT_SECURE_NO_WARNINGS /wd4503 /nologo" )
endif(WIN32)



if( CGAL_FOUND AND Qt5_FOUND AND OPENGL_FOUND )
    include(${CGAL_USE_FILE})
    #include(${QT_USE_FILE})

    SET(MY_BUILD_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build)  #Variable build path
    SET(MY_SRC_PATH ${CMAKE_CURRENT_SOURCE_DIR})  #Variable src path
    INCLUDE_DIRECTORIES(${MY_BUILD_PATH})
    INCLUDE_DIRECTORIES(${MY_SRC_PATH})


    # UI files (Qt Designer files)
    qt5_wrap_ui( DT_UI_FILES delaunay.ui)

    # qrc files (resources files, that contain icons, at least)
    qt5_add_resources( DT_RESOURCE_FILES delaunay.qrc )

    # The executable itself.
    add_executable( delaunay ${SRCS} ${DT_UI_FILES} ${DT_RESOURCE_FILES} )
    qt5_use_modules(delaunay Xml Script OpenGL)

    # Link with Qt libraries
    target_link_libraries( delaunay ${QT_LIBRARIES} )
    
    # Link with Glut and OpenGL
    target_link_libraries( delaunay ${OPENGL_LIBRARY} ${GLUT_LIBRARY}  )

    # Link with CGAL
    target_link_libraries( delaunay ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()
  message(STATUS "NOTICE: This program requires CGAL and Qt5, and will not be compiled.")

endif()
