#include <igl/readOFF.h>
//#include <igl/readPLY.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/implicit_primitives.cpp>
#include <igl/avg_edge_length.h>
#include <igl/writeOFF.h>
#include <iostream>
#include <cmath>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

void signed_distance_sphere(Eigen::MatrixXd& GV, Eigen::VectorXd& S)
{
  // Fill the distance to a unit sphere centered at origin here
  int rows = GV.rows();
  S.resize(rows);
  for (int i = 0; i < rows; ++i) {
    double x = GV(i,0);
    double y = GV(i,1);
    double z = GV(i,2);
    S(i) = 0.1 - std::sqrt(x*x + y*y + z*z);
  }
}

//look carefully textbook and PDF
void grid_spheres(Eigen::MatrixXd& GV, Eigen::VectorXd& S)
{
    int rows = GV.rows();
    S.resize(rows);
    for (int i = 0; i < rows; ++i) {
      double x = GV(i,0);
      double y = GV(i,1);
      double z = GV(i,2);
      //center of sphere 
      double xc = V(0,0);
      double yc = V(0,1);
      double zc = V(0,2);
      
      S(i) = 0.005 - std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
      for(int j = 1 ; j < V.rows() ; j++ ){
	      xc = V(j,0);
	      yc = V(j,1);
	      zc = V(j,2);
      
	      S(i) = std::max( S(i),  0.005 - std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc)));
      }
  }
}

int main( int argc, char *argv[] )
{
  igl::readOFF("/Users/HRMS/Desktop/Gt/data_for_GT/bunny_botsch_simplified_500F.off", V, F);
  double p[4];
  double avg;
  const int s = 128;
  double x_max, y_max, z_max;
  double x_min, y_min, z_min;
  //const Eigen::RowVector3i res(65,65,65);
  //Eigen::MatrixXd GV(res(0)*res(1)*res(2),3);

  x_max = y_max = z_max = -99999.0;
  x_min = y_min = z_min =  99999.0;
  
  for(int i = 0 ; i < V.rows() ; i++ ){
     if( x_max < V(i,0) ) x_max = V(i,0);
	 if( y_max < V(i,1) ) y_max = V(i,1);
	 if( z_max < V(i,2) ) z_max = V(i,2);
	 if( x_min > V(i,0) ) x_min = V(i,0);
	 if( y_min > V(i,1) ) y_min = V(i,1);
	 if( z_min > V(i,2) ) z_min = V(i,2);
  }

  std::cout<<x_max<<", "<<y_max<<", "<<z_max<<std::endl;
  std::cout<<x_min<<", "<<y_min<<", "<<z_min<<std::endl;
  //return 0;

  //calcurate a diagonal
  const double min = -0.2*sqrt(std::pow(x_max-x_min,2.0)+std::pow(y_max-y_min,2.0)+std::pow(z_max-z_min,2.0));
  const double max =  0.2*sqrt(std::pow(x_max-x_min,2.0)+std::pow(y_max-y_min,2.0)+std::pow(z_max-z_min,2.0));


  //enlarged
  //std::cout<<min<<std::endl;
  const Eigen::RowVector3d Vmin(x_min+min,y_min+min,z_min+min);
  const Eigen::RowVector3d Vmax(x_max+max,y_max+max,z_max+max);

  const double h = (Vmax-Vmin).maxCoeff()/(double)s;
  const Eigen::RowVector3i res = (s*((Vmax-Vmin)/(Vmax-Vmin).maxCoeff())).cast<int>();
  
  Eigen::MatrixXd GV(res(0)*res(1)*res(2),3);
 
  //Eigen::VectorXd S;
  //S.resize(GV.rows());
  //p[3]=0.2*igl::avg_edge_length(V,F);
  //for(int i = 0 ; i < GV.rows() ; i++ ){
  //   p[0]=GV(i,0);
  //   p[1]=GV(i,1);
  //   p[2]=GV(i,2);
  //   S(i)=sphere(2.0,2.0,2.0,p);
  //}

  //calculate grid subdivision by x/128 
  for(int zi = 0;zi<res(2);zi++)
  {
    //cartesian look test p.17
    const auto lerp = [&](const int di, const int d)->double
      {return Vmin(d)+(double)di/(double)(res(d)-1)*(Vmax(d)-Vmin(d));};
    const double z = lerp(zi,2);
    for(int yi = 0;yi<res(1);yi++)
    {
      const double y = lerp(yi,1);
      for(int xi = 0;xi<res(0);xi++)
      {
        const double x = lerp(xi,0);
        GV.row(xi+res(0)*(yi + res(1)*zi)) = Eigen::RowVector3d(x,y,z);
      }
    }
  }

  
  Eigen::VectorXd S;
  S.resize(GV.rows());
  grid_spheres(GV,S);
 
  Eigen::MatrixXd SV;
  Eigen::MatrixXi SF; 
  igl::copyleft::marching_cubes(S,GV,res(0),res(1),res(2),SV,SF);

  //v2(2) = V(F(0,1),2);

  //std::cout << V.rows() << std::endl;
  //std::cout << V.cols() << std::endl;

  //igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F);
  //viewer.data().set_face_based(true);
  //viewer.launch();
  bool r = igl::writeOFF("bunny_spheret.off", SV, SF);
  assert(r);

  return 0;
}
