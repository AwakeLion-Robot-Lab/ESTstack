/**
 * \file se2_localization_eskfom.cpp
 *
 *  Created on: Mar 24, 2026
 *     \author: siyiovo
 *
 *  ---------------------------------------------------------
 *  This file is:
 *  (c) 2026 ESTstack Contributors
 *
 *  Licensed under the Apache License, Version 2.0
 *  ---------------------------------------------------------
 *
 *  ---------------------------------------------------------
 *  Demonstration example:
 *
 *  2D Robot localization based on fixed beacons using ESKFOM.
 *
 *  This example is based on the SE2 localization example from
 *  the manif library (se2_localization.cpp) but uses ESTstack's
 *  ESKFOM implementation.
 *
 *  Reference: https://github.com/artivis/manif/blob/devel/examples/se2_localization.cpp
 *  ---------------------------------------------------------
 *
 *  We consider a robot in the plane surrounded by a small
 *  number of punctual landmarks or _beacons_.
 *  The robot receives control actions in the form of axial
 *  and angular velocities, and is able to measure the location
 *  of the beacons w.r.t its own reference frame.
 *
 *  The robot pose X is in SE(2) and the beacon positions b_k in R^2,
 *
 *          | cos th  -sin th   x |
 *      X = | sin th   cos th   y |  // position and orientation
 *          |   0        0      1 |
 *
 *      b_k = (bx_k, by_k)           // lmk coordinates in world frame
 *
 *  The control signal u is a twist in se(2) comprising longitudinal
 *  velocity v and angular velocity w, with no lateral velocity
 *  component, integrated over the sampling time dt.
 *
 *      u = (v*dt, 0, w*dt)
 *
 *  The control is corrupted by additive Gaussian noise u_noise,
 *  with covariance
 *
 *    Q = diagonal(sigma_v^2, sigma_s^2, sigma_w^2).
 *
 *  This noise accounts for possible lateral slippage u_s
 *  through a non-zero value of sigma_s,
 *
 *  At the arrival of a control u, the robot pose is updated
 *  with X <-- X * Exp(u) = X ⊞ u.
 *
 *  Landmark measurements are of the range and bearing type,
 *  though they are put in Cartesian form for simplicity.
 *  Their noise n is zero mean Gaussian, and is specified
 *  with a covariances matrix R.
 *  We notice the rigid motion action y = h(X,b) = X^-1 * b
 *
 *      y_k = (brx_k, bry_k)       // lmk coordinates in robot frame
 *
 *  We consider the beacons b_k situated at known positions.
 *  We define the pose to estimate as X in SE(2).
 *  The estimation error dx and its covariance P are expressed
 *  in the tangent space at X.
 *
 *  The motion and measurement models are
 *
 *    X_(t+1) = f(X_t, u) = X_t * Exp(u)     // motion equation
 *    y_k     = h(X, b_k) = X^-1 * b_k        // measurement equation
 *
 */

#include "eststack/solution/eskfom.hpp"
#include "eststack/model/motion/se2.hpp"
#include "eststack/model/measurement/landmark_se2.hpp"

#include <vector>
#include <iostream>
#include <iomanip>

using std::cout;
using std::endl;

using namespace Eigen;

typedef Array<double, 2, 1> Array2d;
typedef Array<double, 3, 1> Array3d;

int main()
{
    std::srand((unsigned int) time(0));

    // START CONFIGURATION
    //
    //
    const int NUMBER_OF_LMKS_TO_MEASURE = 3;

    // Define the robot pose element and its covariance
    manif::SE2d X_simulation, X_unfiltered;
    Matrix3d    P0;

    X_simulation.setIdentity();
    X_unfiltered.setIdentity();
    P0.setZero();

    // Create ESKFOM filter with options matching manif example (no Joseph form, no reset Jacobian)
    eststack::solution::ESKFOM<manif::SE2d>::Options eskfom_options;
    eskfom_options.use_joseph_form = false;      // Use standard form as in manif example
    eskfom_options.apply_reset_jacobian = false; // No reset Jacobian as in manif example

    manif::SE2d X_initial;
    X_initial.setIdentity();
    eststack::solution::ESKFOM<manif::SE2d> eskf(X_initial, P0, eskfom_options);

    // Define a control vector and its noise and covariance
    manif::SE2Tangentd  u_simu, u_est, u_unfilt;
    Vector3d            u_nom, u_noisy, u_noise;
    Array3d             u_sigmas;
    Matrix3d            U;

    u_nom    << 0.1, 0.0, 0.05;
    u_sigmas << 0.1, 0.1, 0.1;
    U        = (u_sigmas * u_sigmas).matrix().asDiagonal();

    // Define three landmarks in R^2
    Eigen::Vector2d b0, b1, b2, b;
    b0 << 2.0,  0.0;
    b1 << 2.0,  1.0;
    b2 << 2.0, -1.0;
    std::vector<Eigen::Vector2d> landmarks;
    landmarks.push_back(b0);
    landmarks.push_back(b1);
    landmarks.push_back(b2);

    // Define the beacon's measurements
    Vector2d                y, y_noise;
    Array2d                 y_sigmas;
    Matrix2d                R;
    std::vector<Vector2d>   measurements(landmarks.size());

    y_sigmas << 0.01, 0.01;
    R        = (y_sigmas * y_sigmas).matrix().asDiagonal();

    // Create motion and measurement models
    eststack::model::MotionSE2 motion_model;

    //
    //
    // CONFIGURATION DONE

    // DEBUG
    cout << std::fixed   << std::setprecision(3) << std::showpos << endl;
    cout << "X STATE     :    X      Y    THETA" << endl;
    cout << "----------------------------------" << endl;
    cout << "X initial   : " << X_simulation.log().coeffs().transpose() << endl;
    cout << "----------------------------------" << endl;
    // END DEBUG

    // START TEMPORAL LOOP
    //
    //

    // Make 10 steps. Measure up to three landmarks each time.
    for (int t = 0; t < 10; t++)
    {
        //// I. Simulation ###############################################################################

        /// simulate noise
        u_noise = u_sigmas * Array3d::Random();             // control noise
        u_noisy = u_nom + u_noise;                          // noisy control

        u_simu   = u_nom;
        u_est    = u_noisy;
        u_unfilt = u_noisy;

        /// first we move - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        X_simulation = X_simulation.rplus(u_simu);          // X * exp(u)

        /// then we measure all landmarks - - - - - - - - - - - - - - - - - - - -
        for (std::size_t i = 0; i < landmarks.size(); i++)
        {
            b = landmarks[i];                               // lmk coordinates in world frame

            /// simulate noise
            y_noise = y_sigmas * Array2d::Random();         // measurement noise

            y = X_simulation.inverse().act(b);              // landmark measurement, before adding noise
            y = y + y_noise;                                // landmark measurement, noisy
            measurements[i] = y;                            // store for the estimator just below
        }

        //// II. Estimation using ESKFOM ###############################################################

        /// First we move - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        eskf.predict(motion_model, u_est, U);

        /// Then we correct using the measurements of each lmk - - - - - - - - -
        for (int i = 0; i < NUMBER_OF_LMKS_TO_MEASURE; i++)
        {
            // landmark
            b = landmarks[i];                               // lmk coordinates in world frame

            // measurement
            y = measurements[i];                            // lmk measurement, noisy

            // Create measurement model for this landmark
            eststack::model::LandmarkSE2 meas_model(b);

            // Update step
            eskf.update(meas_model, y, R);
        }

        //// III. Unfiltered ##############################################################################

        // move also an unfiltered version for comparison purposes
        X_unfiltered = X_unfiltered.rplus(u_unfilt);

        //// IV. Results ##############################################################################

        // Get estimated state from filter
        auto X_estimated = eskf.getState();

        // DEBUG
        cout << "X simulated : " << X_simulation.log().coeffs().transpose() << endl;
        cout << "X estimated : " << X_estimated.log().coeffs().transpose() << endl;
        cout << "X unfilterd : " << X_unfiltered.log().coeffs().transpose() << endl;
        cout << "----------------------------------" << endl;
        // END DEBUG
    }

    //
    //
    // END OF TEMPORAL LOOP. DONE.

    return 0;
}
