#include <iostream>
#include <vector>
#include <valarray>
#include <memory>
#include <assert.h>

//General Polynomial Spline container
class Spline;
std::shared_ptr<Spline> splinePtr;
class Spline{
  int degree;//coefficients in polynomial segments
  int dimension;//number of individual splines
  double collocation_interval;//time interval between collocation points
  double t0;//initial time
  //coefficient_array[collocation_index][coefficient_index][coordinate_index]
  std::vector< std::vector< std::valarray<double> > > coefficient_array;
public: 
  Spline(const std::vector<std::vector<std::valarray<double>>> _coeff_array=NULL, 
         const double& _collocation_interval, 
         const double& _t0, 
         const double& _dimension, 
         const double& _degree) : 
         collocation_interval(_collocation_interval),
         t0(_t0),
         dimension(_dimension), 
         degree(_degree){
    
    if(not (_coeff_array==NULL)){coefficient_array=_coeff_array;}
  }

  //Copy tail into the back of this Spline
  void concatenate(const splinePtr& tail){
    coefficient_array.push_back(tail->coefficient_array);
  }
  //Push a collocation point into the spline
  void push(const std::vector< std::valarray<double> >& knot){coefficient_array.push_back(knot);}
  
  //Evaluate spline at time t
  std::valarray<double> at(const double& t){
    int index = std::min(coefficient_array.size(),std::max(0,std::floor((t-t0)/collocation_interval)));
    double time = (t-t0)-collocation_interval*((double)index);
    std::valarray<double> eval(0.0,degree);
    for(int i=0;i<degree;i++){
      eval+=coefficient_array[index][i]*pow(t,i);
    }
    return eval;
  }
  // Allocates memory in coefficient_array for "size" collocation points
  void reserve(const int& size){coefficient_array.reserve(size);}
  
};

//Vector field with numerical integration
class DynamicalSystem{
public:
  const double max_time_step;
  int sim_counter=0;
  
  DynamicalSystem(double _max_time_step) : max_time_step(_max_time_step){}
  
  virtual void flow(vctr& dx, const vctr& x,const vctr& u)=0;
  
  virtual double getLipschitzConstant()=0;
  
  //Explicit integration step from x(t1) to x_plus(t2) with input u(t). Assigns traj to spline extending from x
  virtual void step(splinePtr& traj, const std::valarray<double>& x, const splinePtr& u, const double& t1, const double t2)=0;
  
  void sim(splinePtr& solution, double t0, double tf, const std::valarray<double>& x0, const splinePtr& u){
    assert(tf>t0);
    
    //compute minimum number of steps to satisfy dt<dt_max
    double num_steps=ceil((tf-t0)/max_time_step);
    double integration_step=(tf-t0)/num_steps;
    //resize solution
    solution->reserve(num_steps);
    //set initial state and time
    std::valarray<double> state=x0;
    double time=t0;
    splinePtr traj_segment=nullptr;
    //integrate
    for(int i=0;i<num_steps;i++){
      //Use numerical integration scheme to compute a spline extending from state with input u([t,t+integration_step])
      step(traj_segment,state,u,time,time+integration_step);
      //add traj_segment to solution
      solution->concatenate(traj_segment);
      time+=integration_step;
      state=traj_segment->at(time);
    }
    sim_counter++;
    return;
  }
  
};

class SymplecticEuler : public DynamicalSystem{
  std::valarray<double> x1,x2,f0,f1,f2;
  double h;
public:
  //Symplectic or "modified" Euler integration with cubic interpolation between integration steps
  SymplecticEuler(double _max_time_step) : DynamicalSystem(_max_time_step){}
  void step(splinePtr& segment, const std::valarray<double>& x0, const splinePtr& u, const double t1, const double t2) override {  
    
    assert(t1<t2 && "Integration step must be positive in SymplecticEuler");
    
    //Symplectic Euler step
    h=t2-t1;
    flow(f0,x0,u->at(t1));
    x1=x0+0.5*h*f0;
    flow(f1,x1,u->at(t1+0.5*h));
    x2=x0+h*f1;
    flow(f2,x2,u->at(t2));
    
    //Cubic interpolation between x0 and x2 with x'(t1)=f(x0,u(t0)) and x'(t2)=f(x2,u(t2))
    std::vector< std::valarray<double> > cubic;
    cubic.push_back(x0);//t^0 term
    cubic.push_back(f0);//t^1 term
    cubic.push_back((-2*f0+3*f1-f2)/h);//t^2 term
    cubic.push_back((f0-2*f1+f2)/(h*h));//t^3 term
    std::vector<std::vector< std::valarray<double>>> knot_point({cubic});
    segment = new Spline(knot_point,t2-t1,t1,x0.size(),4);
  }
};

class DoubleIntegrator : public SymplecticEuler{
  DoubleIntegrator(double _max_time_step) : SymplecticEuler(_max_time_step){}
  void flow(/TODO)
};

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    return 0;
}
