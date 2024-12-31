cmake_minimum_required(VERSION 3.10)

#Nombre del proyecto
project(project)

# Versión de C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

################ MAIN 
add_library(eigen_lib INTERFACE)

include_directories(/home/rinese/project_ws/lib/project_lib
                        /home/rinese/iris_ws/lib/mosek/10.2/tools/platform/linux64x86/h
                        /home/rinese/project_ws/lib/eigen
                        /usr/include/libqhullcpp
                        /usr/include                        
                        )

add_executable(main 
                        /home/rinese/project_ws/src/main.cpp
                        /home/rinese/project_ws/src/Polyhedron_class.cpp)

target_link_libraries(main PRIVATE
                        /home/rinese/project_ws/lib/mosek/10.2/tools/platform/linux64x86/bin/libmosek64.so.10.2
                        /home/rinese/project_ws/lib/mosek/10.2/tools/platform/linux64x86/bin/libfusion64.so.10.2
                        eigen_lib
                        qhullcpp
                        /home/rinese/project_ws/lib/matplot/lib/libmatplot.a
                        /home/rinese/project_ws/lib/matplot/lib/Matplot++/libnodesoup.a
                        jpeg
                        png
                        z
                        pthread
                        tiff
                        )

################ Inscribed_ellipsoid

#Definimos la biblioteca como PRIVATE si solo queremos que la use este ejecutable

#Especificamos la ruta a los fichero cabecera
include_directories(/home/rinese/project_ws/lib/mosek/10.2/tools/platform/linux64x86/h)

#Directorio con el codigo fuente
file(GLOB SRC_SOURCES ${CMAKE_SOURCE_DIR}/src/inscribed_ellipsoid_n.cpp)

#Crear ejecutable del proyecto
add_executable(inscribed_ellipsoid_n ${SRC_SOURCES})

target_link_libraries(inscribed_ellipsoid_n PRIVATE
        /home/rinese/project_ws/lib/mosek/10.2/tools/platform/linux64x86/bin/libmosek64.so.10.2
        /home/rinese/project_ws/lib/mosek/10.2/tools/platform/linux64x86/bin/libfusion64.so.10.2)


############### Closest_Point_debug

include_directories(/usr/include/libqhullcpp)

#Directorio con el codigo fuente: creamos la variable SRC_SOURCES
file(GLOB SRC_SOURCES ${CMAKE_SOURCE_DIR}/src/closest_point_debug.cpp)

add_executable(closest_point_debug ${SRC_SOURCES})

#Vinculamos la libreria al ejecutable
target_link_libraries(closest_point_debug PRIVATE
        eigen_lib
        qhullcpp
        /home/rinese/project_ws/lib/mosek/10.2/tools/platform/linux64x86/bin/libmosek64.so.10.2
        /home/rinese/project_ws/lib/mosek/10.2/tools/platform/linux64x86/bin/libfusion64.so.10.2
        )

############### Pruebas Matplot

include_directories(/home/rinese/project_ws/lib/matplot/include
			/usr/include/libqhullcpp	
			)

#Directorio con el codigo fuente: creamos la variable SRC_SOURCES
file(GLOB SRC_SOURCES ${CMAKE_SOURCE_DIR}/src/plot_test.cpp)

add_executable(plot_test ${SRC_SOURCES})

#Vinculamos la libreria al ejecutable
target_link_libraries(plot_test PRIVATE
	/home/rinese/project_ws/lib/matplot/lib/libmatplot.a
	/home/rinese/project_ws/lib/matplot/lib/Matplot++/libnodesoup.a
	jpeg
	png
	z
	pthread
	tiff
	eigen_lib
	qhullcpp
        )

############### Matplot

#file(GLOB_RECURSE SOURCES
#	/home/rinese/project_ws/lib/project_lib/GCS_Data.h
#	/home/rinese/project_ws/src/GCS_Data.cpp
#)

include_directories(/home/rinese/project_ws/lib/matplot/include
                        /usr/include/libqhullcpp
			/home/rinese/project_ws/lib/project_lib
			/home/rinese/project_ws/lib/clipper/CPP/Clipper2Lib/include
                        )

#Directorio con el codigo fuente: creamos la variable SRC_SOURCES
file(GLOB SRC_SOURCES ${CMAKE_SOURCE_DIR}/src/plot_graph.cpp)

add_executable(plot_graph ${SRC_SOURCES})

#Vinculamos la libreria al ejecutable
target_link_libraries(plot_graph PRIVATE
        /home/rinese/project_ws/lib/matplot/lib/libmatplot.a
        /home/rinese/project_ws/lib/matplot/lib/Matplot++/libnodesoup.a
        jpeg
        png
        z
        pthread
        tiff
        eigen_lib
        qhullcpp
	/home/rinese/project_ws/lib/clipper/CPP/build/libClipper2.a
	#	/home/rinese/project_ws/src/GCS_Data.cpp
        )

