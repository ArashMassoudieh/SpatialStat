VTKBUILDPATH = /home/arash/Projects/VTK-build
VTKHEADERPATH = /home/arash/Projects/VTK

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++14

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += Include
INCLUDEPATH += ../Utilities
DEFINES += Use_Armadillo

#Adding VTK Libraries
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkChartsCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonColor-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonComputationalGeometry-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonDataModel-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonExecutionModel-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonMath-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonMisc-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonSystem-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkCommonTransforms-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkDICOMParser-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkexpat-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersAMR-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersExtraction-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersFlowPaths-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeneral-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeneric-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersGeometry-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersHybrid-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersHyperTree-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersImaging-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersModeling-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersParallel-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersParallelImaging-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersPoints-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersProgrammable-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSelection-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSMP-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersSources-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersStatistics-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersTexture-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersTopology-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkFiltersVerdict-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkfreetype-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkGeovisCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkgl2ps-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkglew-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkhdf5-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkhdf5_hl-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingColor-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingFourier-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingGeneral-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingHybrid-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingMath-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingMorphological-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingSources-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingStatistics-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkImagingStencil-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInfovisCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInfovisLayout-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionImage-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionStyle-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkInteractionWidgets-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOAMR-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOEnSight-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOExodus-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOGeometry-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOImage-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOImport-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOInfovis-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOLegacy-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOLSDyna-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOMINC-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOMovie-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIONetCDF-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOParallel-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOParallelXML-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOPLY-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOSQL-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOTecplotTable-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOVideo-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOXML-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkIOXMLParser-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkjpeg-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkjsoncpp-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtklibharu-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtklibxml2-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtklz4-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkmetaio-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkParallelCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkpng-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingAnnotation-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingContext2D-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingFreeType-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingGL2PSOpenGL2-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingImage-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingLabel-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingLOD-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingOpenGL2-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingVolume-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkRenderingVolumeOpenGL2-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtksqlite-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtksys-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtktiff-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkverdict-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsContext2D-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsCore-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkViewsInfovis-9.0
LIBS += -L$$VTKBUILDPATH/lib/ -lvtkzlib-9.0
LIBS += -L"/usr/local/lib/ -lsuperlu.so"

#VTK Include files
INCLUDEPATH +=$${VTKHEADERPATH}/Common/Core
INCLUDEPATH +=$${VTKBUILDPATH}/Common/Core
INCLUDEPATH +=$${VTKBUILDPATH}/Common/Color
INCLUDEPATH +=$${VTKBUILDPATH}/Common/DataModel
INCLUDEPATH +=$${VTKBUILDPATH}/Utilities/KWIML
INCLUDEPATH +=$${VTKHEADERPATH}/Utilities/KWIML
INCLUDEPATH +=$${VTKHEADERPATH}/Rendering/Core
INCLUDEPATH +=$${VTKBUILDPATH}/Rendering/Core
INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Core
INCLUDEPATH +=$${VTKHEADERPATH}/Charts/Core
INCLUDEPATH +=$${VTKBUILDPATH}/Charts/Core
INCLUDEPATH +=$${VTKBUILDPATH}/Filters/General
INCLUDEPATH +=$${VTKBUILDPATH}/Rendering/Context2D
INCLUDEPATH +=$${VTKHEADERPATH}/Rendering/Context2D
INCLUDEPATH +=$${VTKHEADERPATH}/Common/DataModel
INCLUDEPATH +=$${VTKHEADERPATH}/Common/Math
INCLUDEPATH +=$${VTKHEADERPATH}/Views/Context2D
INCLUDEPATH +=$${VTKBUILDPATH}/Views/Context2D
INCLUDEPATH +=$${VTKBUILDPATH}/Views/Core
INCLUDEPATH +=$${VTKBUILDPATH}/Interaction/Widgets
INCLUDEPATH +=$${VTKHEADERPATH}/Views/Core
INCLUDEPATH +=$${VTKHEADERPATH}/Interaction/Style
INCLUDEPATH +=$${VTKBUILDPATH}/Interaction/Style
INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Modeling
INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Modeling
INCLUDEPATH +=$${VTKHEADERPATH}/Common/ExecutionModel
INCLUDEPATH +=$${VTKBUILDPATH}/Common/ExecutionModel
INCLUDEPATH +=$${VTKHEADERPATH}/Interaction/Widgets/
INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Core/
INCLUDEPATH +=$${VTKHEADERPATH}/Common/Misc/
INCLUDEPATH +=$${VTKBUILDPATH}/Common/Misc
INCLUDEPATH +=$${VTKHEADERPATH}/IO/XML/
INCLUDEPATH +=$${VTKBUILDPATH}/IO/XML
INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Sources
INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Sources
INCLUDEPATH +=$${VTKHEADERPATH}/Filters/General
INCLUDEPATH +=$${VTKHEADERPATH}/IO/Image
INCLUDEPATH +=$${VTKBUILDPATH}/IO/Image
INCLUDEPATH +=$${VTKHEADERPATH}/Imaging/Core
INCLUDEPATH +=$${VTKBUILDPATH}/Imaging/Core

CONFIG(debug, debug|release) {
    message(Building in debug mode)
    #QMAKE_CXXFLAGS+= -fopenmp
    #QMAKE_LFLAGS +=  -fopenmp
    ! macx: LIBS += -lgomp -lpthread -lgsl
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
    ! macx: LIBS += -lgomp -lpthread -lgsl -larmadillo
    macx: LIBS += -lpthread
}


SOURCES += \
    ../Utilities/Concentrations.cpp \
    ../Utilities/Distribution.cpp \
    ../Utilities/Matrix.cpp \
    ../Utilities/Matrix_arma.cpp \
    ../Utilities/Matrix_arma_sp.cpp \
    ../Utilities/NormalDist.cpp \
    ../Utilities/QuickSort.cpp \
    ../Utilities/Utilities.cpp \
    ../Utilities/Vector.cpp \
    ../Utilities/Vector_arma.cpp \
    Src/command.cpp \
    Src/environment.cpp \
    Src/grid.cpp \
    Src/interface.cpp \
    Src/script.cpp \
    main.cpp \
    mainwindow.cpp

HEADERS += \
    ../Utilities/BTC.h \
    ../Utilities/BTCSet.h \
    ../Utilities/Concentrations.h \
    ../Utilities/Distribution.h \
    ../Utilities/Matrix.h \
    ../Utilities/Matrix_arma.h \
    ../Utilities/Matrix_arma_sp.h \
    ../Utilities/NormalDist.h \
    ../Utilities/QuickSort.h \
    ../Utilities/Utilities.h \
    ../Utilities/Vector.h \
    ../Utilities/Vector_arma.h \
    Include/Structs.h \
    Include/command.h \
    Include/environment.h \
    Include/grid.h \
    Include/interface.h \
    Include/script.h \
    Include/vtk.h \
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

message($$VTKBUILDPATH)
message($${VTKBUILDPATH})
message($$LIBS)


win32:CONFIG(release, debug|release): LIBS += -L$$VTKBUILDPATH/lib/release/ -lvtkChartsCore-9.0
else:win32:CONFIG(debug, debug|release): LIBS += -L$$VTKBUILDPATH/lib/debug/ -lvtkChartsCore-9.0
else:unix: LIBS += -L$$VTKBUILDPATH/lib/ -lvtkChartsCore-9.0

INCLUDEPATH += $$VTKBUILDPATH
DEPENDPATH += $$VTKBUILDPATH
