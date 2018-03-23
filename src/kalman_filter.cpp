#include "kalman_filter.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - (H_ * x_);
	CalculcateValues(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
 
  if(px == 0 && py == 0)
	{
			throw std::logic_error("x and y both are zero, Can not calculate Extended Kalman Filter.");
	}
	
	float rho = sqrt(px*px + py*py);
  float phy = atan2(py, px);	
	float rho_dot = (px*vx + py*vy) / rho;
	
  VectorXd h = VectorXd(3);
  h << rho, phy, rho_dot;
  VectorXd y = z - h;
	
	// Normalize angles
  y[1] = atan2(sin(y[1]), cos(y[1]));
// 	if( y[1] > M_PI )
//     y[1] -= 2. * M_PI;
//   else if( y[1] < -M_PI )
//     y[1] += 2. * M_PI;
	
	CalculcateValues(y);
}

void KalmanFilter::CalculcateValues(const VectorXd &y)
{
	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
	MatrixXd K = PHt * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
