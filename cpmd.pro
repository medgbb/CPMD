TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -llapack -larmadillo
LIBS += -L$$TOP_PWD/hflib/lib -lhartree-fock
LIBS += -lboost_filesystem -lboost_system -lboost_mpi -lboost_serialization
LIBS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile)


INCLUDEPATH += $$TOP_PWD/hflib/include


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



#Copy infiles
copydata.commands = $(COPY_DIR) $$TOP_PWD/hflib/infiles $$OUT_PWD
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata

SOURCES += src/cpmd.cpp

TARGET = cpmd
