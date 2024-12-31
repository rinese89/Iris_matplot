#include <cmath>
#include <matplot/matplot.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include<Eigen/Dense>
#include <limits>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullVertex.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullFacet.h>
#include <libqhullcpp/QhullHyperplane.h>

#include <clipper2/clipper.h>
#include "clipper2/clipper.core.h"
#include "clipper2/clipper.engine.h"
#include "clipper2/clipper.offset.h"
#include "clipper2/clipper.minkowski.h"
#include "clipper2/clipper.rectclip.h"

#include <vector>
#include <cstdlib>

using namespace matplot;
using namespace Clipper2Lib;

std::vector<Eigen::Vector2d> centroids;
std::vector<std::vector<Eigen::Vector2d>> vertices_poligonos;

   std::vector<Eigen::Vector2d> calculateIntersectionPoints(const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
    std::vector<Eigen::Vector2d> points;
    int num_constraints = A.rows();
    //new comment......

    for (int i = 0; i < num_constraints; i++) {
        for (int j = 0; j < num_constraints; j++) {
            Eigen::Matrix2d A_sub;
            A_sub.row(0) = A.row(i);
            A_sub.row(1) = A.row(j);

            Eigen::Vector2d b_sub;
            b_sub << b[i], b[j];

            // Resolver el sistema lineal
            if (A_sub.determinant() != 0) { // Solo resolver si el sistema es resoluble
                Eigen::Vector2d x = A_sub.colPivHouseholderQr().solve(b_sub);

                std::cout << "Punto de interseccion: " << std::endl;
                std::cout << x.x() << "," << x.y() << std::endl;

                // Verificar si satisface todas las restricciones
                //if (((A * x).array() <= b.array()).all()) {
                //    points.push_back(x);
                //    for (const auto& vertex : points) {
                //        std::cout << "Point:" << std::endl;
                //        std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
                //    }
                //    
                //}

                if (((A * x).array() - b.array() <= 0.1 ).all()) {
                    points.push_back(x);
                    for (const auto& vertex : points) {
                        std::cout << "Point:" << std::endl;
                        std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
                    }
                    
                }
            }
        }
    }
                    std::cout << "Point Validation:" << std::endl;
                   for (const auto& vertex : points) {
                        std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
                   }
    return points;
}

// Función para calcular el envolvente convexo usando Qhull
std::vector<Eigen::Vector2d> computeConvexHull(const std::vector<Eigen::Vector2d>& points) {
    
    orgQhull::Qhull qhull;
    // Convertir puntos a un formato plano (array 1D)
    std::vector<double> pointArray;
    for (const auto& p : points) {
        pointArray.push_back(p.x());
        pointArray.push_back(p.y());
    }

    std::cout << "Conversion en array plano" << std::endl;
    
    // Construir el envolvente convexo
    qhull.runQhull("convex_hull", 2, points.size(), pointArray.data(), "Qt"); // Qt: rápido y exacto
    
    std::cout << "Construccion envolvente convexo" << std::endl;

    // Extraer vértices del envolvente convexo
    std::vector<Eigen::Vector2d> hullVertices;
    for (const auto& vertex : qhull.vertexList()) {
        hullVertices.emplace_back(vertex.point()[0], vertex.point()[1]);
    }

    return hullVertices;
}

std::vector<Eigen::Vector2d> store_centroids(const Eigen::Vector2d& centroid){
    std::vector<double> centroid_x;
    std::vector<double> centroid_y;

    centroid_x.push_back(centroid.x());
    centroid_y.push_back(centroid.y());

    centroids.push_back(centroid);


    return centroids;
}

Eigen::Vector2d calculateCentroid(const std::vector<Eigen::Vector2d>& points) {
    Eigen::Vector2d centroid(0, 0);
    for (const auto& point : points) {
        centroid += point;
    }
    centroid /= points.size();
    std::cout << "Centroide: " << std::endl;
    std::cout << centroid << std::endl;
    //double cen_x = centroid[0];
    //double cen_y = centroid[1]; 
    //std::vector<double> centroid_x;
    //std::vector<double> centroid_y;
//
    //centroid_x[0] = cen_x;
    //centroid_y[0] = cen_y;    

    //matplot::scatter(centroid_x, centroid_y, 50);
    //matplot::hold(matplot::on);
    return centroid;
}

std::pair<Eigen::Vector2d, std::vector<Eigen::Vector2d>> sortVerticesByPolarAngle(const std::vector<Eigen::Vector2d>& points){
// Función para ordenar los puntos por ángulo polar
//std::vector<Eigen::Vector2d> sortVerticesByPolarAngle(const std::vector<Eigen::Vector2d>& points) {
    
    Eigen::Vector2d centroid = calculateCentroid(points);

    // Copiar los puntos para ordenarlos
    std::vector<Eigen::Vector2d> sorted_points = points;

    // Ordenar por ángulo polar respecto al centroide
    std::sort(sorted_points.begin(), sorted_points.end(), [&centroid](const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
        double angle_a = atan2(a.y() - centroid.y(), a.x() - centroid.x());
        double angle_b = atan2(b.y() - centroid.y(), b.x() - centroid.x());
        return angle_a < angle_b;
    });

    return make_pair(centroid, sorted_points);
}

std::vector<Eigen::Vector2d> plot_convexset(const Eigen::MatrixXd& A, const Eigen::VectorXd& b){
    
    std::vector<Eigen::Vector2d> points = calculateIntersectionPoints(A, b);
    
    // Calcular el envolvente convexo
    std::vector<Eigen::Vector2d> convexHull = computeConvexHull(points);

    //std::vector<Eigen::Vector2d> sorted_points = sortVerticesByPolarAngle(convexHull);
    std::pair<Eigen::Vector2d, std::vector<Eigen::Vector2d>> vertex_centroid = sortVerticesByPolarAngle(convexHull);

    Eigen::Vector2d centroid = vertex_centroid.first;
    std::vector<Eigen::Vector2d> sorted_points = vertex_centroid.second;

    std::vector<Eigen::Vector2d> centroids = store_centroids(centroid);

    std::vector<double> x_vertices;
    std::vector<double> y_vertices;

    vertices_poligonos.push_back(sorted_points);
    
    // Mostrar los puntos ordenados
    std::cout << "Vértices del polígono en orden:" << std::endl;
    for (const auto& vertex : sorted_points) {
        x_vertices.push_back(vertex.x());
        y_vertices.push_back(vertex.y());
        std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
    }

    // Cerrar el polígono volviendo al primer vértice
    x_vertices.push_back(x_vertices[0]);
    y_vertices.push_back(y_vertices[0]);

    // Graficar los puntos
    scatter(x_vertices, y_vertices, 5); // Tamaño del punto: 50

    // Graficar las líneas
    plot(x_vertices, y_vertices, "r-"); // Líneas rojas
    hold(on);

    // Configurar el gráfico
    title("Polígono generado por restricciones");
    xlabel("x");
    ylabel("y");
    grid(true);

    // Mostrar el gráfico
    show();

    return centroids;
}

void plot_convexset_from_clipper(const std::vector<Eigen::Vector2d>& points){
        
    // Calcular el envolvente convexo
    //std::vector<Eigen::Vector2d> convexHull = computeConvexHull(points);

    //std::vector<Eigen::Vector2d> sorted_points = sortVerticesByPolarAngle(convexHull);
    //std::pair<Eigen::Vector2d, std::vector<Eigen::Vector2d>> vertex_centroid = sortVerticesByPolarAngle(points);

    //Eigen::Vector2d centroid = vertex_centroid.first;
    //std::vector<Eigen::Vector2d> sorted_points = vertex_centroid.second;

    //std::vector<Eigen::Vector2d> centroids = store_centroids(centroid);

    std::vector<double> x_vertices;
    std::vector<double> y_vertices;

    vertices_poligonos.push_back(points);
    
    // Mostrar los puntos ordenados
    std::cout << "Vértices del polígono en orden:" << std::endl;
    for (const auto& vertex : points) {
        x_vertices.push_back(vertex.x());
        y_vertices.push_back(vertex.y());
        std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
    }

    // Cerrar el polígono volviendo al primer vértice
    x_vertices.push_back(x_vertices[0]);
    y_vertices.push_back(y_vertices[0]);

    // Graficar los puntos
    fill(x_vertices, y_vertices); // Tamaño del punto: 50

    // Graficar las líneas
    //plot(x_vertices, y_vertices, "g-"); // Líneas rojas
    hold(on);

    // Configurar el gráfico
    title("Polígono generado por restricciones");
    xlabel("x");
    ylabel("y");
    grid(true);

    // Mostrar el gráfico
    show();

}

void plot_graph(std::vector<Eigen::Vector2d> centroids_vector){

    std::vector<double> x_centroids;
    std::vector<double> y_centroids;

    for (const auto& points : centroids) {
       x_centroids.push_back(points.x());
       y_centroids.push_back(points.y());
    }
    
    
    scatter(x_centroids, y_centroids, 10);
    plot(x_centroids, y_centroids, "r-");

    show();

}

//Clipper2Lib::Path64 convert2ClipperType(const std::vector<Eigen::Vector2d>& vertices) {
//    Clipper2Lib::Path64 path;
//    for (const auto& vertex : vertices) {
//        path.emplace_back(static_cast<int64_t>(vertex.x() * 1000),
//                          static_cast<int64_t>(vertex.y() * 1000));
//    }
//    return path;
//}

PathD convert2ClipperTypeD(const std::vector<Eigen::Vector2d>& vertices) {
    PathD path;
    for (const auto& vertex : vertices) {
        path.emplace_back(vertex.x(), vertex.y());
    }
    return path;
}

// Creamos la union de todos los poligonos que hemos generado con las restricciones
// Leemos de uno en uno los poligonos (vertices_poligonos contiene los poligonos definidos por sus vertices)
// Los converimos a formato PathsD
void create_clipper_polygons(){
    PathsD clipperPoly;
    PathsD clipperPolygons;
    for (const auto& poly : vertices_poligonos) {
        clipperPoly.push_back(convert2ClipperTypeD(poly));
        clipperPolygons = Union(clipperPolygons, clipperPoly, FillRule::NonZero);
    }

    for (int i = 0; i < clipperPolygons.size(); ++i) {
    std::cout << "ClipperPolygons: " << i + 1 << ":" << std::endl;
    for (const auto& vertex : clipperPolygons[i]) {
        std::cout << "(" << vertex << ")" << std::endl;
    }
    std::cout << "---" << std::endl; // Separador entre polígonos
    }

    PathsD boundingBox;

    std::vector<Eigen::Vector2d> bounding_vector = {
        {-10.0, -10.0},
        {10.0, -10.0},
        {10.0, 10.0},
        {-10.0, 10.0}
    };
   
    boundingBox.push_back(convert2ClipperTypeD(bounding_vector));

    PathsD uncoveredRegions; // Contenedor para la diferencia
    uncoveredRegions = Difference(boundingBox, clipperPolygons, FillRule::NonZero);

    //To plot clipper surface: pasamos los tipo PathsD a tipo std::vector<Eigen::Vector2d>
    std::vector<std::vector<Eigen::Vector2d>> poligonos_from_clip;
    
    Eigen::Vector2d point_vertex;

    std::cout << "Regiones no cubiertas (espacios vacíos):" << std::endl;
    for (const auto& region : uncoveredRegions) {
        std::vector<Eigen::Vector2d> vertices_clip;
        std::cout << "Polígono:" << std::endl;
        for (const auto& point : region) {
            point_vertex[0]=point.x;
            point_vertex[1]=point.y;
            vertices_clip.push_back(point_vertex);
            std::cout << "(" << point.x << ", " << point.y << ")" << std::endl;

        }
        poligonos_from_clip.push_back(vertices_clip);
        std::cout << "---" << std::endl;
    }

    //Imprimir la conversion de tipo
    for(const auto& poligon : poligonos_from_clip)
    {
        plot_convexset_from_clipper(poligon);
    for (const auto& vertex : poligon) {
        std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
    }
        std::cout << "----------" << std::endl;
    }
}

//Paths64<int64_t> convertToPaths64(const Paths<int64_t>& paths) {
//    Paths64<int64_t> paths64;
//    for (const auto& path : paths) {
//        Path64<int64_t> path64;
//        for (const auto& point : path) {
//            path64.push_back(Point64(point.X, point.Y));
//        }
//        paths64.push_back(path64);
//    }
//    return paths64;
//}

int main() {
    
    Eigen::Matrix<double, 7, 2> A0({ {-0.187845,-0.982199},{0.989463,-0.144787},{-0.944799,0.327649},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 7> b0({ 3.134441,4.368166,4.472645,10.000000,10.000000,10.000000,10.000000 });
    
    Eigen::Matrix<double, 6, 2> A1({ {-0.447214,-0.894427},{-0.989463,0.144787},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 6> b1({ -9.838699,-4.368166,10.000000,10.000000,10.000000,10.000000 });
    
    Eigen::Matrix<double, 7, 2> A2({ {0.975684,0.219180},{-0.164399,-0.986394},{0.944799,-0.327649},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 7> b2({ -3.782521,3.123581,-4.472645,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 7, 2> A3({ {-0.947879,0.318630},{-1.000000,-0.000000},{0.447214,0.894427},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 7> b3({ -6.945773,-7.000000,9.838699,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 7, 2> A4({ {0.876956,-0.480570},{-0.393919,0.919145},{0.187845,0.982199},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 7> b4({ 3.760376,-2.888742,-3.134441,10.000000,10.000000,10.000000,10.000000 });
	
	Eigen::Matrix<double, 7, 2> A5({ {-0.057434,-0.998349},{-0.975684,-0.219180},{0.944799,-0.327649},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 7> b5({ -5.760362,3.782521,-4.472645,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 8, 2> A6({ {0.987924,0.154942},{-0.934022,-0.357216},{-0.876956,0.480570},{0.187845,0.982199},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 8> b6({ 3.855021,0.989684,-3.760376,-3.134441,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 8, 2> A7({ {-0.707107,-0.707107},{0.057434,0.998349},{-0.975684,-0.219180},{0.944799,-0.327649},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 8> b7({ -1.414214,5.760362,3.782521,-4.472645,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 9, 2> A8({ {-0.903372,-0.428857},{0.987991,-0.154513},{0.707107,-0.707107},{-0.989463,0.144787},{0.554700,0.832050},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 9> b8({ -2.359359,7.066563,7.071068,-4.368166,6.101702,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 7, 2> A9({ {0.000000,1.000000},{-0.876956,0.480570},{0.934022,0.357216},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 7> b9({ -8.000000,-3.760376,-0.989684,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 8, 2> A10({ {0.903372,0.428857},{1.000000,-0.000000},{-0.987924,-0.154942},{-0.989463,0.144787},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 8> b10({ 2.359359,5.000000,-3.855021,-4.368166,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 8, 2> A11({ {1.000000,0.000000},{-0.987924,-0.154942},{-0.918198,0.396121},{0.000000,1.000000},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 8> b11({ 7.000000,-3.855021,-7.363839,-7.000000,10.000000,10.000000,10.000000,10.000000 });
	
    Eigen::Matrix<double, 7, 2> A12({ {-0.707107,0.707107},{-0.707107,-0.707107},{1.000000,-0.000000},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 7> b12({ -7.071068,-1.414214,7.000000,10.000000,10.000000,10.000000,10.000000 });
	
	Eigen::Matrix<double, 7, 2> A13({ {0.928477,-0.371391},{0.393919,-0.919145},{0.164399,0.986394},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 7> b13({ -5.756555,2.888742,-3.123581,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 11, 2> A14({ {-0.164399,-0.986394},{0.187845,0.982199},{-0.944799,0.327649},{0.707107,0.707107},{0.707107,-0.707107},{0.903372,0.428857},{1.000000,0.000000},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 11> b14({ 3.123581,-3.134441,4.472645,1.414214,7.071068,2.359359,7.000000,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 8, 2> A15({ {0.947879,-0.318630},{0.554700,0.832050},{-0.987991,0.154513},{-0.707107,-0.707107},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 8> b15({ 6.945773,6.101702,-7.066563,-1.414214,10.000000,10.000000,10.000000,10.000000 });

	Eigen::Matrix<double, 8, 2> A16({ {0.447214,0.894427},{0.970143,-0.242536},{-0.989463,0.144787},{-0.187845,-0.982199},{-1.000000,0.000000},{1.000000,0.000000},{0.000000,-1.000000},{0.000000,1.000000} });
	Eigen::Vector<double, 8> b16({ 9.838699,3.880570,-4.368166,3.134441,10.000000,10.000000,10.000000,10.000000 });

    plot_convexset(A0,b0);
    plot_convexset(A1,b1);
    plot_convexset(A2,b2);
    plot_convexset(A3,b3);
    plot_convexset(A4,b4);
    plot_convexset(A5,b5);
    plot_convexset(A6,b6);
    plot_convexset(A7,b7);
    plot_convexset(A8,b8);
    plot_convexset(A9,b9);
    plot_convexset(A10,b10);
    plot_convexset(A11,b11);
    plot_convexset(A12,b12);
    plot_convexset(A13,b13);
    plot_convexset(A14,b14);
    plot_convexset(A15,b15);

    //Eigen::Vector2d centroid0 = plot_convexset(A0,b0);
    //Eigen::Vector2d centroid1 = plot_convexset(A1,b1);
    //Eigen::Vector2d centroid2 = plot_convexset(A2,b2);
    //Eigen::Vector2d centroid3 = plot_convexset(A3,b3);
    //Eigen::Vector2d centroid4 = plot_convexset(A4,b4);
    //Eigen::Vector2d centroid5 = plot_convexset(A5,b5);
    //Eigen::Vector2d centroid6 = plot_convexset(A6,b6);
    //Eigen::Vector2d centroid7 = plot_convexset(A7,b7);
    //Eigen::Vector2d centroid8 = plot_convexset(A8,b8);
    //Eigen::Vector2d centroid9 = plot_convexset(A9,b9);
    //Eigen::Vector2d centroid10 = plot_convexset(A10,b10);
    //Eigen::Vector2d centroid11 = plot_convexset(A11,b11);
    //Eigen::Vector2d centroid12 = plot_convexset(A12,b12);
    //Eigen::Vector2d centroid13 = plot_convexset(A13,b13);
    //Eigen::Vector2d centroid14 = plot_convexset(A14,b14);
    //Eigen::Vector2d centroid15 = plot_convexset(A15,b15);
    //Eigen::Vector2d centroid16 = plot_convexset(A16,b16);

    for (const auto& centroid : centroids) {
        std::cout << "(" << centroid.x() << ", " << centroid.y() << ")" << std::endl;
    }

    for (int i = 0; i < vertices_poligonos.size(); ++i) {
        std::cout << "Polígono " << i + 1 << ":" << std::endl;
        for (const auto& vertex : vertices_poligonos[i]) {
            std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
        }
        std::cout << "---" << std::endl; // Separador entre polígonos
    }

    create_clipper_polygons();
    //Clipper2Lib::Paths64 clipperPolygons;
    //for (const auto& poly : vertices_poligonos) {
    //    Clipper2Lib::Path64 clip_polygon = convert2ClipperType(poly);
    //    clipperPolygons.push_back(clip_polygon);
    //}

    
    
   // clipperPolygons.push_back(convert2ClipperTypeD(vertices_poligonos[0]));

    
    

   // Clipper2Lib::Paths64<int64_t> clipperPolygons64 =  convertToPaths64(clipperPolygons);

   //Clipper2Lib::Path64 boundingBox = {
   //    {static_cast<int64_t>(-10000), static_cast<int64_t>(-10000)},
   //    {static_cast<int64_t>(10000), static_cast<int64_t>(-10000)},
   //    {static_cast<int64_t>(10000), static_cast<int64_t>(10000)},
   //    {static_cast<int64_t>(-10000), static_cast<int64_t>(10000)}};

   // hold(off);
   
 
        

    //Unimos todos los poligonos
   // Clipper2Lib::Clipper64 clipperUnion;
   // clipperUnion.AddPaths(clipperPolygons, PathType::Subject, false); 

   // ClipperD clipperUnion;
   // clipperUnion.AddSubject(clipperPolygons);
   // PathsD unionResult, unionResult1;
   // bool union_result_flag = clipperUnion.Execute(ClipType::Union, FillRule::NonZero, unionResult);
//
   // unionResult1 = Union(clipperPolygons, FillRule::NonZero);

    ///////////////////
  //for (int i = 0; i < unionResult1.size(); ++i) {
  //      std::cout << "UNION: " << i + 1 << ":" << std::endl;
  //      for (const auto& vertex : unionResult1[i]) {
  //          std::cout << "(" << vertex << ")" << std::endl;
  //      }
  //      std::cout << "---" << std::endl; // Separador entre polígonos
  //  }
    
   // hold(off);


    ////////////////////




    //ClipperD clipperDifference;
    //clipperDifference.AddSubject({boundingBox});   // Agregar el boundingBox como sujeto principal
    //clipperDifference.AddClip(clipperPolygons);
//
    

    

    //Clipper2Lib::Paths<double> unionResult;
    //clipperUnion.Execute(Union, unionResult);
//
   // // Calcular los espacios vacíos
   // Clipper clipperDifference;
   // clipperDifference.AddPath(boundingBox, ptSubject, true); // Área delimitadora
   // clipperDifference.AddPaths(unionResult, ptClip, true);  // Polígonos combinados
//
   // Paths uncoveredRegions;
   // clipperDifference.Execute(ctDifference, uncoveredRegions);
//
   // // Mostrar los polígonos de las regiones no cubiertas
   // std::cout << "Regiones no cubiertas (espacios vacíos):" << std::endl;
   // for (const auto& region : uncoveredRegions) {
   //     std::cout << "Polígono:" << std::endl;
   //     for (const auto& point : region) {
   //         std::cout << "(" << point.X / 1000.0 << ", " << point.Y / 1000.0 << ")" << std::endl;
   //     }
   //     std::cout << "---" << std::endl;
   // }


    //plot_graph(centroids);
    
    //hold(off);
    //scatter(centroids, centroid.y(), 200)
    //plot(centroid.x(), centroid.y(),  "r-");
    

    return 0;
}
