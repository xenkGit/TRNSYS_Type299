#pragma once
#include <Eigen/Dense>
namespace loadAgg
{
    Eigen::VectorXd time_ClaessonJaved(double timestep, double t_simEnd, int cells_per_level);
    double compute_deltaT(const Eigen::VectorXd& dg, const Eigen::RowVectorXd& q_b);
    Eigen::VectorXd initialize_dg(const Eigen::VectorXd& g_d, double lambda_ground);
    Eigen::MatrixXd initialize_A(const Eigen::VectorXd& width);
    void next_time_step(Eigen::RowVectorXd& q_b, const Eigen::MatrixXd& A);
    void set_current_load(Eigen::RowVectorXd& q_b, double q_b_new);
}