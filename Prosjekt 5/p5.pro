TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    DiffSolv.cpp \
    lib.cpp

HEADERS += \
    DiffSolv.h \
    lib.h


LIBS += -llapack -lblas -larmadillo
