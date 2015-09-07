TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp \
    planet.cpp

LIBS += -armadillo

HEADERS += \
    system.h \
    planet.h
