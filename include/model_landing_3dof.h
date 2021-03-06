#pragma once

#include "constants.hpp"
#include "EcosWrapper.hpp"
#include <Eigen/Dense>
using namespace Eigen;


class model_landing_3dof {


    const double TWR = 2.0;
    const double g = 9.81;
    const double rG = 12.0; // Radius of gyration
    const double rTB = 20.0; // Lever arm, distance: engines to CG

    const double max_gimbal_angle = 0.5;

    double rx_init = 40;
    double ry_init = 100;
    double vx_init = -30;
    double vy_init = 0;
    double theta_init = 0;
    double dtheta_init = 0;

public:

    static constexpr size_t n_states = 6;
    static constexpr size_t n_inputs = 2;
    static string get_name(){ return "model_landing_3dof"; }

    using StateVector   = Eigen::Matrix<double, n_states,        1>;
    using ControlVector = Eigen::Matrix<double, n_inputs,        1>;
    using StateMatrix   = Eigen::Matrix<double, n_states, n_states>;
    using ControlMatrix = Eigen::Matrix<double, n_states, n_inputs>;

    double total_time_guess() { return 8; }

    void initialize(Eigen::Matrix<double, n_states, K> &X, Eigen::Matrix<double, n_inputs, K> &U);

    StateVector                ode(const StateVector &x, const ControlVector &u);
    StateMatrix     state_jacobian(const StateVector &x, const ControlVector &u);
    ControlMatrix control_jacobian(const StateVector &x, const ControlVector &u);

    void add_application_constraints(optimization_problem::SecondOrderConeProgram &socp);

    StateVector get_random_state();
    ControlVector get_random_input();

};