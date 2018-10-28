#include <igl/readOFF.h>
//#include <igl/readPLY.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/implicit_primitives.cpp>
#include <igl/avg_edge_length.h>
#include <iostream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

double average(double, int);

int main( int argc, char *argv[] )
{
  igl::readOFF("/Users/HRMS/Desktop/Gt/data_for_GT/bunny_botsch_simplified_500F.off", V, F);
  //igl::readPLY("/Users/HRMS/Desktop/Gt/data_for_GT/bunny_botsch.ply", V, F);
  //double sum = 0.0;
  //Eigen::Vector3d v1 = Eigen::VectorXd::Zero(3);
  //Eigen::Vector3d v2 = Eigen::VectorXd::Zero(3);
  double p[4];
  double avg;
  double S[V.rows()];
  double x_max, y_max, z_max;
  double x_min, y_min, z_min;

  x_max = y_max = z_max = x_min = y_min = z_min = 0.0;
  
  p[3]=0.2*igl::avg_edge_length(V,F);
  for(int i = 0 ; i < V.rows() ; i++ ){
 	 p[0]=V(i,0);
 	 p[1]=V(i,1);
 	 p[2]=V(i,2);
  	 S[i]=sphere(0.0,0.0,0.0,p);
  }

  for(int i = 0 ; i < V.rows() ; i++ ){
     if( x_max < V(i,0) ) x_max = V(i,0);
	 if( y_max < V(i,1) ) y_max = V(i,1);
	 if( z_max < V(i,2) ) z_max = V(i,2);
	 if( x_min > V(i,0) ) x_min = V(i,0);
	 if( y_min > V(i,1) ) y_min = V(i,1);
	 if( z_min > V(i,2) ) z_min = V(i,2);
  }

  
  igl::copyleft::marching_cubes(S,GV,65,65,65,V,F);

  //v2(2) = V(F(0,1),2);

  //std::cout << V.rows() << std::endl;
  //std::cout << V.cols() << std::endl;

  //igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F);
  //viewer.data().set_face_based(true);
  //viewer.launch();
}
