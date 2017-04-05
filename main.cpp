#include <iostream>
#include <numerical_integration.h>

//Dynamical system using Symplectic Euler integration scheme
class DoubleIntegrator : public SymplecticEuler{
public:
  DoubleIntegrator(double _max_time_step) : SymplecticEuler(_max_time_step){}
  void flow(std::valarray<double>& dx, const std::valarray<double>& x, const std::valarray<double>& u)override{
    assert(x.size()==2 && u.size()==1);
    dx.resize(2);
    dx[0]=x[1];
    dx[1]=u[0];
  }
  double getLipschitzConstant(){return 1.0;}
};

int main(){
  double t0=0.0;//initial time for simulation
  double tf=10.0;//simulation final time
  std::valarray<double> x0({0.0,0.0});//initial state at t0
  double integration_step=1.0;//time step for integration
  DoubleIntegrator system(integration_step);
  int first_order_hold=2;//number of polynomial coefficients for a first order hold u(t)=c0+c1*t
  int signal_dimension=1;//dimension of control input space
  double signal_collocation_interval=5.0;//control is piecwise linear over 5sec intervals
  double knot0=0.0;//u(0)=0.0
  double knot1=1.0;//u(5.0)=1.0
  double knot2=2.0;//u(10.0)=2.0
  //Calculate polynomial coefficients to collocate the individual polynomials in spline
  std::valarray<double> coeff000({knot0});
  std::valarray<double> coeff010({(knot1-knot0)/signal_collocation_interval});
  std::valarray<double> coeff100({knot1});
  std::valarray<double> coeff110({(knot2-knot1)/signal_collocation_interval});
  
  std::vector<std::valarray<double>> segment1({coeff000,coeff010});//u(t)=coeff000+coeff010*t
  std::vector<std::valarray<double>> segment2({coeff100,coeff110});
  
  std::vector<std::vector<std::valarray<double>>> coefficients({segment1,segment2});
  splinePtr signal( new Spline(coefficients,
                               signal_collocation_interval,
                               t0,
                               signal_dimension,
                               first_order_hold));
  
  //Print input signal
  for(double t=0.0;t<10.0;t+=0.1){
    std::cout<< "u(" << t << ")=" << signal->at(t)[0] << std::endl;
  }
  
  //Simulate system and assign solution to trajectory
  splinePtr trajectory;
  system.sim(trajectory,t0,tf,x0,signal);
  for(double t=0.0;t<10.0;t+=0.1){
    std::cout<< "x(" << t << ")=" << "(" <<trajectory->at(t)[0] << "," << trajectory->at(t)[1] << ")" << std::endl;
  }
}
