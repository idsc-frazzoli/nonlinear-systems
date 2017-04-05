# nonlinear-systems
Input signal to trajectory mapping for nonlinear dynamic models

## Input signals
The input signal is a piecewise polynomial function stored in a Spline object 

## Trajectories
The trajectory is also a piecewise polynomial function stored in a Spline object

## Splines
An instance of the Spline class is used for multivariate piecewise polynomial functions. Collocation of the segments or their derivatives is not handled internally in case collocation is not desired (e.g. discontinuous control inputs). A 3D array "coefficient_array" contains the information to retrieve the value of the function at a specific time "t" with the "x.at(t)" function.

## Numerical integration 
Numerical integration is carried out with a modified Euler step which has local error in O(dt^3). The output spline satisfies x(t0)=x_0, x'(t0)=f(x_0,u(t_0)), x(t0+h)=x_0+hf(x_0+0.5hf(x_0,u(t_0)),u(t_0+0.5h)), x'(t0+h)=f(x(t0+h),u(t0+h))
