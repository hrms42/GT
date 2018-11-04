#include <igl/readOFF.h>
//#include <igl/copyleft/cgal/CSGTree.h>
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

static double compute_dot_product(Eigen::Vector3d& A, Eigen::Vector3d& B) {
  return A(0)*B(0) + A(1)*B(1) + A(2)*B(2);
}


static double compute_norm2(Eigen::Vector3d& v) {
  return std::sqrt(compute_dot_product(v,v));
}

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
      
      S(i) = 0.0625 - std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
      for(int j = 1 ; j < V.rows() ; j++ ){
	      xc = V(j,0);
	      yc = V(j,1);
	      zc = V(j,2);
      
	      S(i) = std::max( S(i),  0.0625 - std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc)));
      }
  }
}

double cylinder_capped(double x, double y, double z, Eigen::VectorXd& P){
	Eigen::Vector3d axis_dir(P(0), P(1), P(2));
	Eigen::Vector3d axis_pos(P(3), P(4), P(5));
	double radius = P(6);
	double height = P(7);

	Eigen::Vector3d diff( x-axis_pos(0), y-axis_pos(1), z-axis_pos(2));
	double lamb = compute_dot_product(axis_dir, diff);
	Eigen::Vector3d v( diff(0) - lamb*axis_dir(0), diff(1) - lamb*axis_dir(1), diff(2) - lamb*axis_dir(2));
	double axis_dist = compute_norm2(v);
	double d = axis_dist - radius;
	d = -d;

	double planes = height/2.0 - std::fabs(lamb);
	double dist = std::min(planes, d);

	return dist;
}

void test_gc(Eigen::MatrixXd& GV, Eigen::VectorXd& S)
{
		int rows = GV.rows();
		S.resize(rows);
 		double height = 3;
		double radius = 0.05;
		int j = 0, k = 0;
		Eigen::VectorXd P(8);

		

		for (int i = 0 ; i < rows ; ++i ){
				double x = GV(i,0);
				double y = GV(i,1);
				double z = GV(i,2);

				P << 1,0,0,0,0,0,radius,height;

				S(i) = cylinder_capped(x,y,z,P);
		}
}


void grid_cylinders(Eigen::MatrixXd& GV, Eigen::VectorXd& S)
{
		int rows = GV.rows();
		S.resize(rows);
 		double height = 1;
		double radius = 0.05;
		int j = 0, k = 0;
		Eigen::VectorXd P(8);

		for (int i = 0 ; i < rows ; ++i ){
				double x = GV(i,0);
				double y = GV(i,1);
				double z = GV(i,2);

				height = std::max( std::fabs(V(F(0,1),0)-V(F(0,0),0)), std::max( std::fabs(V(F(0,1),0)-V(F(0,0),0)), std::fabs(V(F(0,1),0)-V(F(0,0),0))));

                P << (V(F(0,1),0)-V(F(0,0),0)), (V(F(0,1),1)-V(F(0,0),1)), (V(F(0,1),2)-V(F(0,0),2)),
				     (V(F(0,1),0)+V(F(0,0),0))/2.0, (V(F(0,1),1)+V(F(0,0),1))/2.0, (V(F(0,1),2)+V(F(0,0),2))/2.0,
				     radius,height; 

				S(i) = cylinder_capped(x,y,z,P);

		       	for( j = 0 ; j < F.rows(); j++ ){
				   for( k = 1 ; k < F.cols(); k++ ){

				 height = std::max( std::fabs(V( F(j,((k+1)%F.cols())), 0) - V(F(j,k),0)), std::max( std::fabs(( V( F(j,((k+1)%F.cols())), 1) - V(F(j,k),1) )), std::fabs(( V( F(j,((k+1)%F.cols())),2) - V(F(j,k),2) ))));
		   
				 P << ( V( F(j,((k+1)%F.cols())), 0) - V(F(j,k),0) ),  ( V( F(j,((k+1)%F.cols())), 1) - V(F(j,k),1) ), ( V( F(j,((k+1)%F.cols())),2) - V(F(j,k),2) ),
				      ( V( F(j,((k+1)%F.cols())), 0) + V(F(j,k),0) ) / 2.0,  ( V( F(j,((k+1)%F.cols())), 1) + V(F(j,k),1) ) / 2.0, ( V( F(j,((k+1)%F.cols())),2) + V(F(j,k),2) ) / 2.0,
					  radius, height;

					 S(i) = std::max( S(i), cylinder_capped(x,y,z,P)) ;
					 //S(i) = cylinder_capped(x,y,z,P);
				   }
				}
		}
}

int main( int argc, char *argv[] )
{
  igl::readOFF("/Users/HRMS/Desktop/Gt/data_for_GT/simpleBox.off", V, F);
  const int s = 64;
  //double x_max, y_max, z_max;
  //double x_min, y_min, z_min;
  Eigen::Vector3d m = V.colwise().minCoeff();
  Eigen::Vector3d M = V.colwise().maxCoeff();

  //const Eigen::RowVector3i res(65,65,65);
  //Eigen::MatrixXd GV(res(0)*res(1)*res(2),3);

  //x_max = y_max = z_max = -99999.0;
  //x_min = y_min = z_min =  99999.0;
  //
  //for(int i = 0 ; i < V.rows() ; i++ ){
  //       if( x_max < V(i,0) ) x_max = V(i,0);
  //   if( y_max < V(i,1) ) y_max = V(i,1);
  //   if( z_max < V(i,2) ) z_max = V(i,2);
  //   if( x_min > V(i,0) ) x_min = V(i,0);
  //   if( y_min > V(i,1) ) y_min = V(i,1);
  //   if( z_min > V(i,2) ) z_min = V(i,2);
  //}

  //std::cout<<x_max<<", "<<y_max<<", "<<z_max<<std::endl;
  //std::cout<<x_min<<", "<<y_min<<", "<<z_min<<std::endl;

  //calcurate a diagonal
  const double min = -0.2*sqrt(std::pow(M(0)-m(0),2.0)+std::pow(M(1)-m(1),2.0)+std::pow(M(2)-m(2),2.0));
  const double max =  0.2*sqrt(std::pow(M(0)-m(0),2.0)+std::pow(M(1)-m(1),2.0)+std::pow(M(2)-m(2),2.0));

  //enlarged
  //std::cout<<min<<std::endl;
  const Eigen::RowVector3d Vmin(m(0)+min,m(1)+min,m(2)+min);
  const Eigen::RowVector3d Vmax(M(0)+max,M(1)+max,M(2)+max);

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

  //Eigen::VectorXd SC;
  //SC.resize(GV.rows());
  //grid_spheres(GV,SC);
 
  Eigen::VectorXd S;
  S.resize(GV.rows());
  grid_cylinders(GV,S);

  //igl::copyleft::cgal::CSGTree<MatrixXi> CSGTree = {{SC,F},{S,F}};
 
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
  bool r = igl::writeOFF("simpleBoxOut.off", SV, SF);
  assert(r);

  return 0;
}
