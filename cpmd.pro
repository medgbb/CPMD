TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += src/cpmd.cpp

LIBS += -llapack -larmadillo
LIBS += -L$$TOP_PWD/hflib/lib -lhartree-fock


#QMAKE_LIBDIR += $$TOP_PWD/library/HF/lib
INCLUDEPATH += $$TOP_PWD/hflib/include


QMAKE_CXXFLAGS += -std=c++0x


# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc


#C++11 features
QMAKE_CXXFLAGS += -std=c++0x -DARMA_USE_CXX11
QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS -g
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS -O3 -DARMA_NO_DEBUG
