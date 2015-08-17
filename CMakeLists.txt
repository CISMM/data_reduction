cmake_minimum_required(VERSION 2.8)

PROJECT(DataReduction)

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_SOURCE_DIR}
	/usr/include/qwt
	/usr/local/qwt-6.0.0/include
	../extern_libs/qwt6/include
	../extern_libs/jama125
	../extern_libs/tnt_126
)

SET(Sources
cluster.cpp
cluster.h
decompress.cpp
dr_data.h
dr_functions.cpp
file_list.cpp
file_list.h
GLWidget.cpp
qt_test_QImage.cpp
video_analysis.cpp
)
SET(Qt_headers
GLWidget.h
)

QT4_WRAP_CPP(MOCSrcs ${Qt_headers})

SOURCE_GROUP("Generated" FILES
${MOCSrcs}
)

ADD_EXECUTABLE(DataReduction ${Sources} ${MOCSrcs})
TARGET_LINK_LIBRARIES( DataReduction
	qwt
	GLU
	gsl
	gslcblas
)
