#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  x_ = F_*x_;
  P_ = F_*P_*F_ + Q_.transpose();
}

void KalmanFilter::Update(const VectorXd &z) {
  //std::cout << "Enter laser update" << std::endl;
  //std::cout << "z:" << std::endl << z << std::endl;
  //std::cout << "H_*x_:" << std::endl << H_*x_ << std::endl;
  VectorXd y = z-H_*x_;
  //std::cout << "y:" << std::endl << y << std::endl;

  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  x_ += K*y;
  P_ = (MatrixXd::Identity(4,4)-K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //std::cout << "Enter radar update" << std::endl;
  Tools tool;
  VectorXd y(3);
  VectorXd h(3);
  MatrixXd Hj = tool.CalculateJacobian(x_);
  //std::cout << "Hj: " << std::endl << Hj << std::endl;
  double c1 = std::sqrt(x_(0)*x_(0) +x_(1)*x_(1));
  h << c1, std::atan2(x_(1),x_(0)), (x_(0)*x_(2)+x_(1)*x_(3))/c1;
  y = z-h;
  MatrixXd S = Hj*P_*Hj.transpose() + R_;
  MatrixXd K = P_*Hj.transpose()*S.inverse();
  x_ += K*y;
  P_ = (MatrixXd::Identity(4,4)-K*Hj)*P_;
}
