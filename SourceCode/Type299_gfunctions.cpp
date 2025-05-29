#include "Type299_gfunctions.h"

#include <numbers>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_interp.h>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

const double pi = std::numbers::pi;

// collection of functions for g-function calculation based on Düber 2022
namespace gFunc
{
    double ierf(double X)
    {
        return X * std::erf(X) - (1. / std::sqrt(pi)) * (1. - std::exp(-X * X));
    }

    double ils(double h, double d)
    {
        return 2 * ierf(h) + 2 * ierf(h + 2 * d) - ierf(2 * h + 2 * d) - ierf(2 * d);
    }

    double g_ILS(double r, double lambda_ground, double a_ground, double t)
    {
        // Compute the argument for gsl_sf_expint_E1
        double argument = r * r / (4.0 * a_ground * t);
        if (argument > 100) // At the beginning ICS is more important
            return (1.0 / (4.0 * pi * lambda_ground)) * (exp(-argument) / argument);
        else
            return (1.0 / (4.0 * pi * lambda_ground)) * gsl_sf_expint_E1(argument);
    }

    // Struct for integration parameters
    struct FLSParams 
    {
        double r, l_BHE, a_ground, t;
    };

    //Function to integrate
    double integrand(double s, void* p)
    {
        auto* params = static_cast<FLSParams*>(p);
        double r = params->r, l_BHE = params->l_BHE, a_ground = params->a_ground, t = params->t;
        double exponent = -r * r * s * s;
        double ils_value = ils(l_BHE * s, 1 * s); // D=1
        double denominator = l_BHE * s * s;
        double result = std::exp(exponent) * ils_value / denominator;
        return result;
    }

    // g_FLSjc function using numerical integration
    double g_FLSjc(double r, double l_BHE, double lambda_ground, double a_ground, double t, gsl_integration_workspace* w)
    {
        double s = 1. / std::sqrt(4. * a_ground * t);
        gsl_function F;
        FLSParams params = {r, l_BHE, a_ground, t};
        F.function = &integrand;
        F.params = &params;
        double result = 0., error = 0.;
        int status = gsl_integration_qagiu(&F, s, 1.e-15, 1.e-7, w->limit, w, &result, &error);
        return (1. / (4. * pi * lambda_ground)) * result;
    }

    struct ICSParams
    {
        double r;        // Radial distance from the origin to the evaluation point (m).
        double r_BHE;     // Radius of the cylindrical source (borehole radius) (m).
        double a_ground; // Thermal diffusivity of the ground (m^2/s).
        double t;        // Time for which the g-function is being evaluated (s).
    };

    double ics_integrand(double theta, void *p)
    {
        auto* params = static_cast<ICSParams*>(p);
        double r = params->r, r_BHE = params->r_BHE, a_ground = params->a_ground, t = params->t;
        double cos_theta = std::cos(theta);
        double argument = (r * r + r_BHE * r_BHE - 2.0 * r * r_BHE * cos_theta) / (4.0 * a_ground * t);
        if (argument > 100) // TODO: betrifft ca. die ersten 2 Monate
            return (1./pi) * (exp(-argument) / argument); // TODO!!!
        else
            return (1./pi) * gsl_sf_expint_E1(argument); //underflow error
    }

    double g_ICS(double r, double r_BHE, double lambda_ground, double a_ground, double t, gsl_integration_workspace* w)
    {
        gsl_function F;
        ICSParams params = {r, r_BHE, a_ground, t};
        F.function = &ics_integrand;
        F.params = &params;
        double result = 0.0, error = 0.0;
        gsl_integration_qags(&F, 0., pi, 1.e-15, 1.e-7, w->limit, w, &result, &error); 
        return (1.0 / (4.0 * pi * lambda_ground)) * result;
    }

    double g_FullScale(double r, double r_BHE, double l_BHE, double lambda_ground, double a_ground, double t, gsl_integration_workspace* w, std::ofstream& gfsfile)
    { // g-function with FLS corrected for cylindrical geometry with ICS and ILS
        // short- and medium-term: ILS approx.= FLS -> ICS, medium and long term: ICS approx.= ILS -> FLS
        double gICS = g_ICS(r, r_BHE, lambda_ground, a_ground, t, w);
        double gFLS = g_FLSjc(r, l_BHE, lambda_ground, a_ground, t, w);
        double gILS = g_ILS(r, lambda_ground, a_ground, t);
        if (gfsfile.is_open()) // Test: löschen oder als Info lassen?
        {
            gfsfile << t << "\t" << gICS << "\t" << gFLS << "\t" << gILS << "\n";
        }
        return gICS + gFLS - gILS;
    }

    Eigen::VectorXd g_Borefield(
        const Eigen::VectorXd& xPos,
        const Eigen::VectorXd& yPos,
        const Eigen::VectorXd& time,
        double r_BHE, double l_BHE, double lambda_ground, double cvol_ground)
    {
        size_t nBHE = xPos.size();
        size_t nTime = time.size();
        Eigen::VectorXd gFunction(nTime); // Initialize borefield-wide g-function
        double a_ground = lambda_ground / cvol_ground; // thermal diffusivity (m2/s)
        // Compute the distance matrix between boreholes
        Eigen::MatrixXd gMatrixTemp(nBHE, nBHE);
        for (size_t i = 0; i < nBHE; ++i)
        {
            for (size_t j = 0; j < nBHE; ++j)
            {
                double r = sqrt(pow(xPos[i] - xPos[j], 2) + pow(yPos[i] - yPos[j], 2));
                gMatrixTemp(i, j) = (r == 0) ? (r_BHE) : r;
            }
        }
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
        if (!w)
        {
            //logError("Failed to allocate GSL workspace");
            return {};
        }       
        std::ofstream gfsfile("gfunctions.txt"); // Test: löschen oder als Info lassen?
        gfsfile << "Time_s" << "\t" << "gICS" << "\t" << "gFLS" << "\t" << "gILS" << "\n";
        std::ofstream gffile("gfunction.txt");
        gffile << "Time_s" << "\t" << "gFunction" << "\n";
        // Compute average g-function for the entire borefield at each time step
        for(size_t t = 0; t < nTime; ++t)
        {
            double g_sum = 0.0;
            for(size_t i = 0; i < nBHE; ++i)
            {
                for(size_t j = 0; j < nBHE; ++j)
                    g_sum += g_FullScale(gMatrixTemp(i, j), r_BHE, l_BHE, lambda_ground, a_ground, time(t), w, gfsfile);
            }
            gFunction[t] = g_sum / nBHE; // Compute the average g-function
            gffile << time(t) << "\t" << gFunction[t] << "\n";
        }
        gfsfile.close();
        gffile.close();
        gsl_integration_workspace_free(w);
        return gFunction;
    }
}