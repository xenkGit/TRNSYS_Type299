#pragma once
#include <Eigen/Dense>
#include <vector>

namespace gFunc
{
    std::vector<double> geomspace(double start, double stop, int num);
    double ierf(double X);
    double ils(double h, double d);
    double g_FLSjc(double x, double y, double l_BHE, double lambda_ground, double cvol_ground, double t);
    double g_ILS(double x, double y, double lambda_ground, double cvol_ground, double t);
    double g_ICS(double x, double y, double lambda_ground, double cvol_ground, double t);
    double g_FullScale(double x, double y, double l_BHE, double lambda_ground, double cvol_ground, double t);
    
    Eigen::VectorXd g_Borefield(
        const Eigen::VectorXd& xPos,
        const Eigen::VectorXd& yPos,
        const Eigen::VectorXd& time,
        double r_BHE, double l_BHE, double lambda_ground, double cvol_ground);

    /*std::vector<std::vector<std::vector<double>>> g_Matrix(
        const Eigen::VectorXd& xPos,
        const Eigen::VectorXd& yPos,
        const std::vector<double>& time,
        double r_BHE, double l_BHE, double lambda_ground, double cvol_ground);
    Eigen::VectorXd spatiallySuperpose_gFunctions(std::vector<std::vector<std::vector<double>>>& gMatrix, int nBHE, int nTime);*/
}