TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    eig_jac.cpp \
    init_two_elec.cpp \
    init_one_elec.cpp

LIBS += -llapack -lblas -armadillo

HEADERS +=
