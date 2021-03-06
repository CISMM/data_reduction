######################################################################
# Automatically generated by qmake (2.01a) Mon May 21 17:26:42 2012
######################################################################

# This project has been built under Linux successfully by Alfred Zhong.
# It does not currently compile under Windows because it finds undefined
# symbols in Qwt and gsl.

# *** This program uses GPL code, so can only be distributed under GPL.
# Required libraries:
#  jama, tnt (include-only libraries found in ../extern_libs).
#  qw, cluster (need to be built and copied, in ../extern_libs_src).
#  The Gnu Scientific Library (gsl), also included in ../extern_libs_src.
#   On Windows/Cygwin:
#       configure --prefix=/cygdrive/c/usr/local
#       make -j
#       make install
#  To get things going on Ubuntu:
#   sudo apt-get install qtcreator libqwt-dev libgsl0-dev libglu-dev

QT += core gui opengl
# this is the place for switching into C++ 11
QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp

TEMPLATE = app
TARGET = 
DEPENDPATH += .
INCLUDEPATH += .
INCLUDEPATH += ./jama125
INCLUDEPATH += ./tnt_126
INCLUDEPATH += /usr/local/qwt-6.1.0/
INCLUDEPATH += /usr/local/qwt-6.1.0/lib/
INCLUDEPATH += /usr/local/qwt-6.1.0-rc3/include
INCLUDEPATH += /usr/local/qwt-6.1.0/include
INCLUDEPATH += ../extern_libs/qwt6/include
INCLUDEPATH += C:/Qwt-6.0.2-svn/include
INCLUDEPATH += /afs/cs.unc.edu/home/cshao/src/dr_run/extern_libs/jama125
INCLUDEPATH += /afs/cs.unc.edu/home/cshao/src/dr_run/extern_libs/tnt_126
INCLUDEPATH += C:/usr/local/include
INCLUDEPATH += ~/build/include

LIBS += -L~/build/libs
# remove rc3 in lqwt will make it compile on panoptes
#LIBS += -L/usr/local/qwt-6.1.0/lib/ -L../extern_libs/qwt6/lib/ -LC:/Qwt-6.0.2-svn/lib/  -lqwt -L/usr/local/qwt-6.1.0-rc3/lib/
LIBS += -L/usr/local/qwt-6.1.0/lib/ -L../extern_libs/qwt6/lib/ -LC:/Qwt-6.0.2-svn/lib/
#LIBS += -LC:/usr/local/lib/ -lgsl -lgslcblas
LIBS += -lGLU

# Input
HEADERS += dr_data.h GLWidget.h tnt.h tnt_array1d.h tnt_array1d_utils.h tnt_array2d.h  \
    tnt_array2d_utils.h tnt_array3d.h tnt_array3d_utils.h tnt_cmat.h tnt_fortran_array1d.h \
    tnt_fortran_array1d_utils.h tnt_fortran_array2d.h tnt_fortran_array2d_utils.h \
    tnt_fortran_array3d.h tnt_fortran_array3d_utils.h tnt_i_refvec.h tnt_math_utils.h \
    tnt_sparse_matrix_csr.h tnt_stopwatch.h tnt_subscript.h tnt_vec.h tnt_version.h \
    jama_cholesky.h jama_eig.h jama_lu.h jama_qr.h jama_svd.h
SOURCES += dr_functions.cpp GLWidget.cpp video_analysis.cpp qt_test_QImage.cpp \ 
    file_list.cpp \
    main_functions.cpp

CONFIG += console
#QT -= gui
