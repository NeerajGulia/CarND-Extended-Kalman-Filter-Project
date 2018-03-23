#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::endl;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        			0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
							0, 0.0009, 0,
							0, 0, 0.09;

	ekf_.F_ = MatrixXd(4, 4);

	// Process noise covariance matrix
	ekf_.Q_ = MatrixXd(4, 4);
	//create a 4D state vector, we don't know yet the values of the x state
	ekf_.x_ = VectorXd(4);
	
	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
						 0, 1, 0, 0,
						 0, 0, 1, 0,
						 0, 0, 0, 1;
	
	//measurement matrix
	H_laser_ << 1, 0, 0, 0,
			  		0, 1, 0, 0;

	//set the acceleration noise components
	noise_ax = 9;
	noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
  	
		previous_timestamp_ = measurement_pack.timestamp_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
			float rho = measurement_pack.raw_measurements_[0];
			float phy = measurement_pack.raw_measurements_[1];
			float rho_dot = measurement_pack.raw_measurements_[2];
			
			float x = rho * cos(phy);
			float y = rho * sin(phy);
			float vx = rho_dot * cos(phy);
			float vy = rho_dot * sin(phy);
			ekf_.x_ << x, y, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

	/*****************************************************************************
   *  Prediction
   ****************************************************************************/
	//compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;
	
  /* Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds. */
	ekf_.F_ << 1, 0, dt, 0,
	          0, 1, 0, dt,
	          0, 0, 1, 0,
	          0, 0, 0, 1;
	
	//2. Set the process covariance matrix Q
	float dt2 = dt * dt;
	float dt3 = dt2 * dt;
	float dt4 = dt3 * dt;
	float dt3_2 = dt3/2.0;
	float dt4_4 = dt4/4.0;
		
	// Update the process noise covariance matrix.
	ekf_.Q_ << dt4_4 * noise_ax, 0, dt3_2 * noise_ax, 0,
						 0, dt4_4 * noise_ay, 0, dt3_2 * noise_ay,
						 dt3_2 * noise_ax, 0, dt2 * noise_ax, 0,
						 0, dt3_2 * noise_ay, 0, dt2 * noise_ay;
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			ekf_.R_ = R_radar_;
			ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
			ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	}
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			ekf_.R_ = R_laser_;
			ekf_.H_ = H_laser_;
			ekf_.Update(measurement_pack.raw_measurements_);
  }
}
