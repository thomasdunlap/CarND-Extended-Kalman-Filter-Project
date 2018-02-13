#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() != ground_truth.size() || estimations.size() == 0) {
    std::cout << "Error: estimation == 0 or vector size not equal to ground_truth vector size." << std::endl;
    return rmse;
  }

  for (unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];

    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  MatrixXd Hj(3, 4);

  //Check for small values of position magnitude to avoid division by zero
  float rho = pow((pow(px,2) + pow(py,2)), 0.5);
  if(rho < 0.0001){
    cout << "Rho too small. Rho now equals 0.0001";
    rho = 0.0001;
  }

  float inv_rho = pow(rho, -1);
  Hj << px*inv_rho, py*inv_rho, 0, 0,
        -py * pow(inv_rho,2), px * pow(inv_rho, 2), 0, 0,
        py * (vx*py - vy*px) * pow(inv_rho, 3), px * (vy*px - vx*py) * pow(inv_rho, 3), px*inv_rho, py*inv_rho;

  return Hj;
}
