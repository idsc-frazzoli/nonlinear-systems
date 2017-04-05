#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <valarray>
#include <memory>
#include <assert.h>

class Spline;
typedef std::shared_ptr<Spline> splinePtr;
typedef std::vector< std::vector< std::valarray<double> > > splineTensor;

class Spline{
  int degree;//coefficients in polynomial segments
  int dimension;//number of individual splines
  double collocation_interval;//time interval between collocation points
  double t0;//initial time
  
  //coefficient_array[time_interval_index][polynomial_coefficient_index][polynomial_coordinate_index]
  splineTensor coefficient_array;
public: 
  Spline(const splineTensor& _coeff_array, 
         const double& _collocation_interval, 
         const double& _t0, 
         const double& _dimension, 
         const double& _degree) : 
         collocation_interval(_collocation_interval),
         t0(_t0),
         dimension(_dimension), 
         degree(_degree),
         coefficient_array(_coeff_array){}
         
         Spline(const double& _collocation_interval, 
                const double& _t0, 
                const double& _dimension, 
                const double& _degree) : 
                collocation_interval(_collocation_interval),
                t0(_t0),
                dimension(_dimension), 
                degree(_degree){coefficient_array = splineTensor();}
                
                //Copy tail into the back of this Spline 
                void concatenate(const splinePtr& tail){
                  coefficient_array.insert(coefficient_array.end(),
                                           tail->coefficient_array.begin(),
                                           tail->coefficient_array.end());
                }
                //Push a collocation point into the spline
                void push(const std::vector< std::valarray<double> >& knot){coefficient_array.push_back(knot);}
                
                //Evaluate spline at time t
                std::valarray<double> at(const double& t){
                  int index = std::min( (int)coefficient_array.size()-1,std::max(0,(int)std::floor((t-t0)/collocation_interval)));
                  double time = (t-t0)-collocation_interval*((double)index);
                  std::valarray<double> eval(0.0,dimension);
                  for(int i=0;i<degree;i++){
                    eval+=coefficient_array[index][i]*pow(time,i);
                  }
                  return eval;
                }
                // Allocates memory in coefficient_array for "size" collocation points
                void reserve(const int& size){
                  assert(size>=0 && "Cannot reserve negative space for Spline coefficients");
                  coefficient_array.reserve(size);}
                  
};
#endif