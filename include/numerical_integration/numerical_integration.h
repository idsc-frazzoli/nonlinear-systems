#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H

#include <spline.h>

//Vector field with numerical integration
class DynamicalSystem{
protected:
  const double max_time_step;
  int sim_counter=0;
public:
  
  DynamicalSystem(double _max_time_step) : max_time_step(_max_time_step){}
  
  //This is the function defining the dynamic model x'(t)=f(x(t),u(t))
  virtual void flow(std::valarray<double>& dx, const std::valarray<double>& x,const std::valarray<double>& u)=0;
  
  //Planner requires access to a Lipschitz constant on the vector field
  virtual double getLipschitzConstant()=0;
  
  //Assigns segment to a polyniomial connecting x(t1)=x0 to x(t2) 
  virtual void step(splinePtr& segment, const std::valarray<double>& x0, const splinePtr& u, const double& t1, const double& t2)=0;
  
  //Calls step repeatedly to create a spline approximating the solution to an ODE 
  void sim(splinePtr& solution, double t0, double tf, const std::valarray<double>& x0, const splinePtr& u){
    assert(tf>t0);
    double num_steps=ceil((tf-t0)/max_time_step);
    double integration_step=(tf-t0)/num_steps;
    if(solution==nullptr){
      solution = splinePtr(new Spline(integration_step,t0,x0.size(),4));
    }
    solution->reserve(num_steps);
    //set initial state and time
    std::valarray<double> state=x0;
    double time=t0;
    splinePtr traj_segment;
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


//An integration scheme which has local error of O(dt^3)  
class SymplecticEuler : public DynamicalSystem{
protected:
  std::valarray<double> x1,x2,f0,f1,f2;
  double h;
public:
  //Symplectic or "modified" Euler integration with cubic interpolation between integration steps
  SymplecticEuler(double _max_time_step) : DynamicalSystem(_max_time_step){}
  //Step returns a spline between collocation points of the integrator
  void step(splinePtr& segment, const std::valarray<double>& x0, const splinePtr& u, const double& t1, const double& t2) override {  
    assert(t1<t2 && "Integration step must be positive in SymplecticEuler");
    
    //Modified Euler step: local error=O((t2-t1)^3)
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
    cubic.push_back((-2.0*f0+3.0*f1-f2)/h);//t^2 term
    cubic.push_back((f0-2.0*f1+f2)/(h*h));//t^3 term
    std::vector<std::vector< std::valarray<double>>> knot_point({cubic});
    segment=splinePtr(new Spline(knot_point,t2-t1,t1,x0.size(),4));
  }
};

#endif