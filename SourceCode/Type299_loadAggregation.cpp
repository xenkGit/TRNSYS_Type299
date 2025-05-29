#include "Type299_loadAggregation.h"

#include <cmath>
#include <vector>
#include <numbers>
#include <tuple>
#include <Eigen/Dense>

#include <TRNSYS.h>

const double pi = std::numbers::pi;

// collection of functions for load aggregation by Claesson & Javed 2012 based on pygfunction (Cimmino 2018)
namespace loadAgg
{
    Eigen::VectorXd time_ClaessonJaved(double timestep, double t_simEnd, int cells_per_level)
    { // Build a time vector of expanding cell width
        std::vector<double> time_values;
        time_values.reserve(100); // TODO
        double t = 0.0;
        int i = 0;

        while (t < t_simEnd)
        {
            i++;
            double v = std::ceil(static_cast<double>(i) / cells_per_level);
            double width = std::pow(2.0, v - 1);
            t += width * timestep;
            time_values.push_back(t);
        }
        Eigen::VectorXd time_vector = Eigen::Map<Eigen::VectorXd>(time_values.data(), time_values.size());
        return time_vector;
    }

    double compute_deltaT(const Eigen::VectorXd& dg, const Eigen::RowVectorXd& q_b)
    {
        return dg.dot(q_b.transpose());
    }

    Eigen::VectorXd initialize_dg(const Eigen::VectorXd& g_d, double lambda_ground)
    { // Create a matrix of thermal response factor increments for later use in temporal superposition.
        int nt_gd = g_d.size(); //number of time values at which the thermal response factors are required
        Eigen::VectorXd dg(nt_gd);

        double scaling_factor = 1.; //1.0 / (2.0 * pi * lambda_ground); already contained in g_d
        dg(0) = g_d(0) * scaling_factor;
        for (int i = 1; i < nt_gd; ++i)
        {
            double test1 = g_d(i);
            double test2 = g_d(i-1);
            dg(i) = (g_d(i) - g_d(i - 1)) * scaling_factor;
        }
        return dg;
    }

    Eigen::MatrixXd initialize_A(const Eigen::VectorXd& width)
    {
        int nt = width.size();
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nt, nt);
        for (int i = 0; i < nt; ++i)
        {
            double test1 = width(i);
            A(i, i) = 1.0 - 1.0 / width(i);
        }
        // Set superdiagonal elements (1 / width[i+1])
        for (int i = 0; i < nt - 1; ++i)
        {
            double test2 = width(i+1);
            A(i, i + 1) = 1.0 / width(i + 1);
        }
        return A;
    }

    void next_time_step(Eigen::RowVectorXd& q_b, const Eigen::MatrixXd& A)
    { // Shifts aggregated loads by one time step
        q_b *= A;
    }

    void set_current_load(Eigen::RowVectorXd& q_b, double q_b_new)
    { // Set the load at the current time step
        q_b(0) = q_b_new;
    }
}