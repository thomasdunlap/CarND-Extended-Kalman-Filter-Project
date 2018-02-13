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
  Hj << 0,0,0,0,
        0,0,0,0,
        0,0,0,0;

  //Check for small values of position magnitude to avoid division by zero
  float rho = pow((pow(px,2) + pow(py,2)), 0.5);
  if( rho < 0.0001){
    cout << "Rho too small. Rho now equals 0.0001";
    rho = 0.0001;
  }

  float inv_rho = pow(rho, -1);
  Hj(0, 0) = px * inv_rho;
  Hj(1, 0) = -py * pow(inv_rho,2);
  Hj(2, 0) = py * (vx*py - vy*px) * pow(inv_rho, 3);
  Hj(0, 1) = py * inv_rho;
  Hj(1, 1) = px * pow(inv_rho, 2);
  Hj(2, 1) = px * (vy*px - vx*py) * pow(inv_rho, 3);
  Hj(2, 2) = Hj(0, 0);
  Hj(2, 3) = Hj(0, 1);

  

/*
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = c1 * c2;

  if (fabs(c1) < 0.0001) {
    std::cout << "Error: CalculateJacobian() division by zero." << std::endl;
    return Hj;
  }

  Hj << px/c2, py/c2, 0, 0,
        -(py/c1), px/c1, 0, 0,
        py*(vx*py - vy*px) /c3, px*(px*vy - py*vx) /c3, px/c2, py/c2;
*/
  return Hj;
}
