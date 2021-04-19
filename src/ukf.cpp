#include "ukf.h"
#include "Eigen/Dense"

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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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
  

  // set state dimension
  n_x_ = 5;

  // set augmented dimension
  n_aug_ = 7; 

  // set Sigma point spreading parameter
  lambda_ = 3 - n_aug_;


  // create augmented mean vector
  x_aug = VectorXd(7);

  // set sigma point matrix aug
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented state covariance
  P_aug = MatrixXd(7, 7);

  
  // set Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; ++i) // 2n+1 weights
  {  
    weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // set sigma point matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);


  // acc to ukf.h line 45
  is_initialized_ = false;

  noise_radar = MatrixXd(3,3);
  noise_radar <<  std_radr_*std_radr_, 0, 0,
                  0, std_radphi_*std_radphi_, 0,
                  0, 0,std_radrd_*std_radrd_;

  noise_lidar = MatrixXd(2,2);
  noise_lidar <<  std_laspx_*std_laspx_, 0,
                  0, std_laspy_*std_laspy_;                  



  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 0.0225, 0,
        0, 0, 0, 0, 0.0225;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_)  
  {
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);
      x_ << x , y, v, 0, 0;
      
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      x_ << x, y, 0., 0, 0;
      
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    return;

  }

  delta_t_ = (meas_package.timestamp_ - time_us_) / 1000000.0;

  time_us_ = meas_package.timestamp_;

  // Predict
  Prediction(delta_t_);

  // Measurement updates
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }

}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */


  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug_.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1) = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    // extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v_ = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    // predicted state values
    double p_x_pred, p_y_pred , v_pred , yaw_pred , yawd_pred ;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        p_x_pred = p_x + v_/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        p_y_pred = p_y + v_/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        p_x_pred = p_x + v_*delta_t*cos(yaw);
        p_y_pred = p_y + v_*delta_t*sin(yaw);
    }

    v_pred = v_;
    yaw_pred = yaw + yawd*delta_t;
    yawd_pred = yawd;

    // add noise
    p_x_pred = p_x_pred + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    p_y_pred = p_y_pred + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_pred = v_pred + nu_a*delta_t;

    yaw_pred = yaw_pred + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_pred = yawd_pred + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,i) = p_x_pred;
    Xsig_pred_(1,i) = p_y_pred;
    Xsig_pred_(2,i) = v_pred;
    Xsig_pred_(3,i) = yaw_pred;
    Xsig_pred_(4,i) = yawd_pred;
  }


  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) // iterate over sigma points
  {  
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) // iterate over sigma points
  {  
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  
  
  n_z_lidar = 2;
  z_lidar = meas_package.raw_measurements_;

  z_sig_lidar = MatrixXd(n_z_lidar, 2*n_aug_+1);

  for(int i = 0; i < 2*n_aug_+1; i++)
  {
     z_sig_lidar(0, i) = Xsig_pred_(0, i);
     z_sig_lidar(1, i) = Xsig_pred_(1, i);
  }

  z_pred_lidar = VectorXd(n_z_lidar);
  z_pred_lidar.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    z_pred_lidar = z_pred_lidar + weights_(i) * z_sig_lidar.col(i);
  } 

  S_lidar =  MatrixXd(n_z_lidar,n_z_lidar);
  S_lidar.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    VectorXd z_diff = z_sig_lidar.col(i) - z_pred_lidar;
    S_lidar = S_lidar + weights_(i) * z_diff * z_diff.transpose();
  }

  S_lidar = S_lidar + noise_lidar;  // add measurement noise covariance matrix


  Tc_lidar = MatrixXd(n_x_,n_z_lidar);
  Tc_lidar.fill(0.0); 

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    VectorXd z_diff = z_sig_lidar.col(i) - z_pred_lidar;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc_lidar = Tc_lidar + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  K_lidar = Tc_lidar * S_lidar.inverse();

  // residual
  VectorXd z_diff = z_lidar - z_pred_lidar;

  // update state mean and covariance matrix
  x_ = x_ + K_lidar * z_diff;
  P_ = P_ - K_lidar*S_lidar*K_lidar.transpose(); 
  

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  
  n_z_radar = 3;
  z_radar = meas_package.raw_measurements_;

  z_sig_radar = MatrixXd(n_z_radar, 2*n_aug_+1);
  z_sig_radar.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    double px_p_= Xsig_pred_(0, i);
    double py_p_ = Xsig_pred_(1, i);
    double vp_ = Xsig_pred_(2, i);
    double yaw_p_ = Xsig_pred_(3, i);
    double yawd_p_ = Xsig_pred_(4, i);
    double vx_p_ = cos(yaw_p_)*vp_;
    double vy_p_ = sin(yaw_p_)*vp_;

    z_sig_radar(0, i) = sqrt(px_p_*px_p_ + py_p_*py_p_);                      
    z_sig_radar(1, i) = atan2(py_p_, px_p_);                              
    z_sig_radar(2, i) = (px_p_*vx_p_ + py_p_*vy_p_)/(sqrt(px_p_*px_p_ + py_p_*py_p_)); 
  }

  z_pred_radar = VectorXd(n_z_radar);
  z_pred_radar.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    z_pred_radar = z_pred_radar + weights_(i) * z_sig_radar.col(i);
  } 

  S_radar =  MatrixXd(n_z_radar,n_z_radar);
  S_radar.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {  
    VectorXd z_diff = z_sig_radar.col(i) - z_pred_radar;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_radar = S_radar + weights_(i) * z_diff * z_diff.transpose();
  }

  S_radar = S_radar + noise_radar;  // add measurement noise covariance matrix


  Tc_radar = MatrixXd(n_x_,n_z_radar);
  Tc_radar.fill(0.0); 

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    VectorXd z_diff = z_sig_radar.col(i) - z_pred_radar;
    // angle normalization
    while(z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc_radar = Tc_radar + weights_(i) * x_diff * z_diff.transpose();
  }


  // Kalman gain K;
  K_radar = Tc_radar * S_radar.inverse();

  // residual
  VectorXd z_diff = z_radar - z_pred_radar;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K_radar * z_diff;
  P_ = P_ - K_radar*S_radar*K_radar.transpose(); 
  
}