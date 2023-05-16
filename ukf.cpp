#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // flag for first measurement
  is_initialized_ = false;

  // initial state vector
  n_x_ = 5;
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  //P_ = MatrixXd(n_x_, n_x_);
  //P_.fill(0.0);
  P_ = MatrixXd::Identity(n_x_, n_x_);
  //std::cout<<P_<<endl;

  // Initial predict sigma pts matrix 
  n_aug_ = 7; 
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1); 
  Xsig_pred_.fill(0.0);
  
  // initial weights of sigma pts 
  lambda_ = 3 - n_x_;
  weights_ = VectorXd(2*n_aug_ + 1); 
  weights_(0) = lambda_/(lambda_ + n_aug_);
  for (int i =1; i <2*n_aug_+1; ++i)
  {
    weights_(i) = 1/(2*lambda_ + 2*n_aug_);
  }


  // Process noise standard deviation longitudinal acceleration in m/s^2.(modified)
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2 (modified)
  std_yawdd_ = 2.0;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_)
  {
      if (meas_package.sensor_type_ == MeasurementPackage::LASER)
      {
        x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1),0,0,0,0;
        P_(0,0) = std_laspx_ * std_laspx_; 
        P_(1,1) = std_laspy_ * std_laspy_; 
       }
      else if (meas_package.sensor_type_== MeasurementPackage::RADAR)
      {
        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        double rho_dot = meas_package.raw_measurements_[2];

        x_ << rho * cos(phi),  rho*sin(phi),0,0,0;
        P_(0,0) = std_radr_ * std_radr_; 
        P_(1,1) = std_radr_ * std_radr_; 
        P_(2,2) = std_radr_ * std_radr_; 
        P_(3,3) = std_radphi_* std_radphi_;
        P_(4,4) = std_radphi_ * std_radphi_; 
    
      }
      time_us_ = meas_package.timestamp_; 
      is_initialized_ = true; 
      return;

  }

  // time elapse 
  const float dt = (meas_package.timestamp_ - time_us_)/1e6;
  //update timestamp
  time_us_ = meas_package.timestamp_; 

  // Predict the next state 
  Prediction(dt); 

  //Update next state
  if(meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if ( meas_package.sensor_type_== MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else 
  {
    std::cout<<" Invalid sensor info: Cannot update measurment"<<std::endl;
  }

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // Create mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_; 
  x_aug(5) = 0; 
  x_aug(6) =0;

  // Create augmented state cov
  MatrixXd P_aug = MatrixXd::Identity(n_aug_, n_aug_);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5)= std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_; 

  // create sigma pts 
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);

  // create sqrt matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma pts 
  Xsig_aug.col(0) = x_aug; 
  for( int i = 1; i< n_aug_; ++i)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_)* L.col(i);
    Xsig_aug.col(i + n_aug_ +1) = x_aug - sqrt(lambda_ + n_aug_)*L.col(i);
  }

  // Predict sigma pts 
  for(int i = 0; i < 2*n_aug_ +1; ++i)
  {
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predic state 
    double px_pred, py_pred ;

    if (fabs(yawd)> 1e-6)
    {
      px_pred = px + v/yawd* (sin(yaw + yawd *delta_t) - sin(yaw));
      py_pred = py + v/yawd * (-1 * cos(yaw + yawd*delta_t)+ cos(yaw)); 
    }
    else 
    {
      px_pred = px + v*cos(yaw) * delta_t; 
      py_pred = py + v*sin(yaw)* delta_t;
    }
     //add noise for px, py
    px_pred += 0.5 * delta_t*delta_t*cos(yaw) * nu_a; 
    py_pred += 0.5 * delta_t*delta_t*sin(yaw) * nu_a;

    double v_pred = v + delta_t * nu_a; 
    double yaw_pred = yaw + yawd*delta_t + 0.5* delta_t*delta_t*nu_yawdd; 
    double yawd_pred = yawd + delta_t*nu_yawdd;

    // update predict sigma pt matrix 

    Xsig_pred_(0,i) = px_pred; 
    Xsig_pred_(1,i) = py_pred; 
    Xsig_pred_(2,i) = v_pred; 
    Xsig_pred_(3,i) = yaw_pred; 
    Xsig_pred_(4,i) = yawd_pred;
  }


  // Predict mean state

  x_.fill(0.0); 
  for (int i = 0; i <2*n_aug_+1; ++i)
  {
    x_ = x_ + weights_(i)*Xsig_pred_.col(i);       
  }

  // Predict cov
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i)
    {
      VectorXd x_diff = Xsig_pred_.col(i) - x_; 
      //angle normalization
      while(x_diff(3) > M_PI) x_diff(3) -= 2.0 * M_PI; 
      while( x_diff(3) < -M_PI) x_diff(3) += 2.0*M_PI;
      P_ = P_ + weights_(i) * x_diff*x_diff.transpose();
    }
   
}


void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  int n_z_ =2; 
  //create matrix for sigma pts in meas package
  MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_+1);

  // Mean predicted meas.
  VectorXd z_pred = VectorXd(n_z_);
  // Meas. cov matrix 
  MatrixXd S = MatrixXd(n_z_, n_z_);
  
  // Lidar meas 
  VectorXd z = meas_package.raw_measurements_;
  // Cross correlation matrix Tc 
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  // Transform sigma poits to meas. domain 
  for (int i =0; i< 2*n_aug_+1; ++i)
  {
    Zsig(0,i) = Xsig_pred_(0,i); 
    Zsig(1,i) = Xsig_pred_(1,i);
  }

  // Calcualte mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred; 
    //angle normalization 
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI; 
    S += weights_(i)*z_diff*z_diff.transpose();

  }

  // add meas. noise 
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R.fill(0.0); 
  R(0,0) = std_laspx_*std_laspx_; 
  R(1,1) = std_laspy_ * std_laspy_;

  S =  S+ R;

  // Calculate cross correlation matrix
  Tc.fill(0.0); 
  for (int i = 0; i< 2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred; 
    //angle normalization 
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI; 

    VectorXd x_diff = Xsig_pred_.col(i) - x_; 
    while (x_diff(3)>M_PI) x_diff(3) -= 2.0*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.0*M_PI;

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  // Calcualte Kalman gain

  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred; 
  //angle normalization 
  while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI; 

  // update state mean and cov matrix 
  x_ += K*z_diff; 
  P_ -= K*S*K.transpose();

} // end of fnc

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  // Intialize
  int n_z_ = 3; 
  MatrixXd Zsig = MatrixXd(n_z_, 2* n_aug_+1); 
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S = MatrixXd(n_z_, n_z_);
  VectorXd z = meas_package.raw_measurements_; 
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  // transform sigma pts to meas. domain 
  for (int i = 0; i< 2*n_aug_+1;++i )
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double vx  = cos(yaw)*v; 
    double vy = sin(yaw)*v;

    Zsig(0,i) = sqrt(px*px + py*py);
    Zsig(1,i) = atan2(py, px);
    Zsig(2,i) = (px*vx + py*vy)/sqrt(px*px + py*py);
  }

  // Calcualte mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred; 
    //angle normalization 
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI; 
    S += weights_(i)*z_diff*z_diff.transpose();
  }

  // add mea noise cov matrix 
  MatrixXd R = MatrixXd(n_z_, n_z_); 
  R.fill(0.0); 
  R(0,0) = std_radr_ * std_radr_; 
  R(1,1) = std_radphi_*std_radphi_; 
  R(2,2) = std_radrd_ * std_radrd_;

  S +=R; 

  // Calculate cross correlation matrix
  Tc.fill(0.0); 
  for (int i = 0; i< 2*n_aug_+1; ++i)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred; 
    // angle normalization 
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI; 

    VectorXd x_diff = Xsig_pred_.col(i) - x_; 
    while (x_diff(3)>M_PI) x_diff(3) -= 2.0*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.0*M_PI;

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  // Calcualte Kalman gain

  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred; 
  // angle normalization 
  while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI; 

  // update state mean and cov matrix 
  x_ += K*z_diff; 
  P_ -= K*S*K.transpose();

}