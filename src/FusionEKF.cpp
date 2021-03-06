#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // Initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);


  // Laser measurement covariance
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // Radar measurement convariance
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // Laser measurement matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.F_ = MatrixXd(4, 4); // Transition matrix
  ekf_.P_ = MatrixXd(4, 4); // Radar measurement matrix
  ekf_.Q_ = MatrixXd(4, 4); // Process covariance matrix
  /*

  */

  // Acceleration noise in x and y planes
  noise_ax = 9;
  noise_ay = 9;

}

/**
* Destructor.
* Currently nothing here.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
/*
    float px = 0;
    float py = 0;
    float vx = 0;
    float vy = 0;
*/
    // first measurement
    cout << "EKF: " << endl;

    // State vector
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 0.5, 0.5; // Can adjust init values.

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      ekf_.x_[0] = rho * cos(phi);
      ekf_.x_[1] = rho * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
/*
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
*/
      ekf_.x_[0] = measurement_pack.raw_measurements_[0];
      ekf_.x_[1] = measurement_pack.raw_measurements_[1];
    }


    previous_timestamp_ = measurement_pack.timestamp_;

    // Initial customized covariance matrix
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 500, 0,
               0, 0, 0, 500;

    // Initial state transition matrix
    ekf_.F_ << 1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;


    is_initialized_ = true;
    cout << "Initalized FusionEKF." << '\n';
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**

     * Update the state transition matrix F according to the new elapsed time.
      - dt is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  // Acceleration noise components


  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  ekf_.Q_ << dt_4/4 * noise_ax, 0, dt_3/2 * noise_ax, 0,
             0, dt_4/4 * noise_ay, 0, dt_3/2 * noise_ay,
             dt_3/2 * noise_ax, 0, dt_2 * noise_ax, 0,
             0, dt_3/2 * noise_ay, 0, dt_2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;

    ekf_.R_ = MatrixXd(3, 3);
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.H_ = H_laser_;

    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
