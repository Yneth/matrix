#-------------------------------------------------
#
# Project created by QtCreator 2015-03-18T22:03:50
#
#-------------------------------------------------

QT       -= gui

TARGET = Matrix
TEMPLATE = lib
CONFIG += staticlib

SOURCES += matrix.cpp\
        vector.cpp\
        augmentedmatrix.cpp\
    main.cpp

HEADERS += matrix.h\
        vector.h\
        augmentedmatrix.h
unix {
    target.path = /usr/lib
    INSTALLS += target
}
