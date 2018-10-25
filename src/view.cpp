#include <igl/readOFF.h>
//#include <igl/readPLY.h>
#include <igl/opengl/glfw/Viewer.h>

#include <iostream>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main( int argc, char *argv[] )
{
  igl::readOFF("/Users/HRMS/Desktop/Gt/data_for_GT/bunny_botsch_simplified_500F.off", V, F);
  //igl::readPLY("/Users/HRMS/Desktop/Gt/data_for_GT/bunny_botsch.ply", V, F);

  std::cout << V.rows() << std::endl;
  std::cout << V.cols() << std::endl;
  std::cout << F.rows() << std::endl;
  std::cout << F.cols() << std::endl;
 
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.data().set_face_based(true);
  viewer.launch();
}
