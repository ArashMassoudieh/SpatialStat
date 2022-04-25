QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += Include
INCLUDEPATH += ../Utility_Classes

CONFIG(debug, debug|release) {
    message(Building in debug mode)
    #QMAKE_CXXFLAGS+= -fopenmp
    #QMAKE_LFLAGS +=  -fopenmp
    ! macx: LIBS += -lgomp -lpthread
    macx: LIBS += -lpthread
    DEFINES += NO_OPENMP DEBUG

} else {
    message(Building in release mode)
    !macx:QMAKE_CXXFLAGS += -fopenmp
    !macx:QMAKE_LFLAGS +=  -fopenmp
    # QMAKE_CFLAGS+=-pg
    # QMAKE_CXXFLAGS+=-pg
    # QMAKE_LFLAGS+=-pg
    macx: DEFINES += NO_OPENMP
    ! macx: LIBS += -lgomp -lpthread
    macx: LIBS += -lpthread
}


SOURCES += \
    ../Utility_Classes/DistributionNUnif.cpp \
    ../Utility_Classes/Matrix.cpp \
    ../Utility_Classes/Matrix_arma.cpp \
    ../Utility_Classes/NormalDist.cpp \
    ../Utility_Classes/QuickSort.cpp \
    ../Utility_Classes/Utilities.cpp \
    ../Utility_Classes/Vector.cpp \
    ../Utility_Classes/Vector_arma.cpp \
    Src/command.cpp \
    Src/grid.cpp \
    Src/interface.cpp \
    main.cpp \
    mainwindow.cpp

HEADERS += \
    ../Utility_Classes/BTC.h \
    ../Utility_Classes/BTCSet.h \
    ../Utility_Classes/DistributionNUnif.h \
    ../Utility_Classes/Matrix.h \
    ../Utility_Classes/Matrix_arma.h \
    ../Utility_Classes/NormalDist.h \
    ../Utility_Classes/QuickSort.h \
    ../Utility_Classes/Utilities.h \
    ../Utility_Classes/Vector.h \
    ../Utility_Classes/Vector_arma.h \
    Include/Structs.h \
    Include/command.h \
    Include/grid.h \
    Include/interface.h \
    mainwindow.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

win32 {

    LAPACK_INCLUDE = $$PWD/include
    #64 bits build
    contains(QMAKE_TARGET.arch, x86_64) {
        #debug
        CONFIG(debug, debug|release) {
            LAPACK_LIB_DIR = $$PWD/libs/lapack-blas_lib_win64/debug
            LIBS +=  -L$${LAPACK_LIB_DIR} -llapack_win64_MTd \
                    -lblas_win64_MTd
        }
        #release
        CONFIG(release, debug|release) {
            LAPACK_LIB_DIR = $$PWD/libs/lapack-blas_lib_win64/release
            LIBS +=  -L$${LAPACK_LIB_DIR} -llapack_win64_MT \
                    -lblas_win64_MT
        }
    }

    INCLUDEPATH += $${LAPACK_INCLUDE}

    DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS

}

linux {
    #sudo apt-get install libblas-dev liblapack-dev
     DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS
     LIBS += -larmadillo -llapack -lblas
}

macx {
    #sudo apt-get install libblas-dev liblapack-dev
     DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS
     LIBS += -llapack -lblas
}
