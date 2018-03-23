#include <iostream>
#include "tools.h"
#include <cmath>
#include <stdexcept>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::pow;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
  rmse << 0,0,0,0;
    
  if(estimations.size() == 0)
  {
      throw std::logic_error("empty estimations.");
  }
  if(estimations.size() != ground_truth.size() )
  {
      throw std::logic_error("Estimations vector length should be equal to ground truth vectors.");
  }
	//accumulate squared residuals
  VectorXd residual;
  for(int i=0; i < estimations.size(); ++i){
        residual = ground_truth[i] - estimations[i];
        residual = residual.array()*residual.array();
        rmse += residual;
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE 
  float c = (px * px) + (py * py);
  float c_sqrt = sqrt(c);
  float c_3by2 = c * c_sqrt; //pow(c, 3.0/2.0);

	//check division by zero
  if(c == 0)
  {
    throw std::logic_error("division by zero");
  }
	//initialize 
	Hj << 0, 0, 0, 0,
	      0, 0, 0, 0,
      	0, 0, 0, 0;

  //compute the Jacobian matrix
  Hj(0,0) = px/c_sqrt;
  Hj(0,1) = py/c_sqrt;
  Hj(1,0) = -py/c;
  Hj(1,1) = px/c;
  Hj(2,0) = py*(vx*py - vy*px)/c_3by2;
  Hj(2,1) = px*(vy*px - vx*py)/c_3by2;
  Hj(2,2) = Hj(0,0);
  Hj(2,3) = Hj(0,1);

	return Hj;
}