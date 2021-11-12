# adaptive-estimation
Multiple model adaptive estimation for state and parameter estimation

In this repository, I applied the multiple model adaptive estimation on a differential drive robot model to track the states as well as the process noise covariance.
This parallel algorithm is applicable for nonlinear systems with changes in the motion model as well as the sensor uncertainties. For more detailed information about this filter and its uderlying mathematics, see "Optimal Estimation for Dynamic Systems", John L. Crassidis, John L. Junkins  

Here is a brief description of the functions in the code:

function xdot=diff_drive(x,u)

This function runs the differential drive motion model with a state vector x and a controller u.


function [t,ym]=generate_measurement(tspan,x_0,u,h,r,q_truth)

This method generates a series of artificial measurements to feed to the system as the observation data.
The parameters are time span tspan, initial condition x_0, controller u, the sensor model matrix h, measurement noise covariance r, and the true value of the process noise covariance denoted by q_truth. The output is a time vector and the measurement ym


function [F,b]=linearize1(x_0,u)

This function linearizes the process model


function xdot=f(t,xx,u,coeff,q)

This method runs a deterministic version of the differential drive robot.
