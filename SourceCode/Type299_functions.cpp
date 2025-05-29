#include <cmath>
#include <vector>
#include <numbers>
#include <tuple>

#include "Type299_functions.h"
#include <TRNSYS.h>

const double pi = std::numbers::pi;

namespace t299
{
    void message(int crntUnit, int crntType, const char errorType[], const char message[])
    { // Pass string message to Trnsys
        int errorCode = -1;
        char type[20];
        char message4Trnsys[400];
        strcpy_s(type,errorType);
        strcpy_s(message4Trnsys,message);
        messages(&errorCode, message4Trnsys, type, &crntUnit, &crntType, (size_t)strlen(message4Trnsys), (size_t)strlen(type));
    }

    double calcNusseltNo(double re, double pr, int crntUnit, int crntType)
    { // Calculate the Nusselt number
        double nu;
        if((3000 <= re) && (re <= pow(10,6))) // Taler & Taler 2017
            {
            if((0.1 <= pr) && (pr <= 1.))
            {
                nu = 0.02155*pow(re,0.8018)*pow(pr,0.7095);
            }
            else if ((1. <= pr) && (pr <= 3.))
            {
                nu = 0.01253*pow(re,0.8413)*pow(pr,0.6179);
            }
            else if ((3. <= pr) && (pr <= 1000.))
            {
                nu = 0.00881*pow(re,0.8991)*pow(pr,0.3911);
            }
            else
            {
                nu = 0.00881*pow(re,0.8991)*pow(pr,0.3911);
                t299::message(crntUnit, crntType, "Warning", "Prandtl number is outside the acceptable range.");
            }
        }
        else if (re < 3000.) // TODO decide on 2300 (laminar, Laferrière) or 3000 (turbulent, Taler²) - used also for stagnation
        {
            nu = 3.66; // Laferrière 2020 bzw. VDI Wärmeatlas: Laminare thermisch und hydrodynamisch ausgebildete Strömung durch Rohre bei konstanter Wandtemperatur
        }
        else
        {
            nu = 0.00881*pow(re,0.8991)*pow(pr,0.3911);
            t299::message(crntUnit, crntType, "Warning", "Reynolds number is outside the acceptable range.");
        }
        
        return nu;
    }

    std::tuple<double,double,double>calcGroutResistances(double res_fp, double lambda_grout, double lambda_ground, double r_BHE, double r_tube, double d_shanks, int n_pipes)
    { // Calculate Rb and Ra using explicit multipole formulas according to Claesson and Javed 2019
        double sigma = (lambda_grout - lambda_ground) / (lambda_grout + lambda_ground);
        double res_b0 = res_fp/4 + 1/(8*pi*lambda_grout) * (std::log(pow(r_BHE,4) / (4*r_tube*pow(d_shanks,3))) + sigma*std::log(pow(r_BHE,8) / (pow(r_BHE,8) - pow(d_shanks,8))));
        
        double beta = 2*pi*lambda_grout*res_fp;
        double b1 = (1-beta) / (1+beta);
        double ppc = pow(r_tube,2) / (4*pow(d_shanks,2));
        double pc = pow(d_shanks,2) / pow(pow(r_BHE,8)-pow(d_shanks,8),1/4);
        double pb = pow(r_BHE,2) / pow(pow(r_BHE,8) - pow(d_shanks,8),1/4);
        
        double res_b1 = res_b0 - 1/(8*pi*lambda_grout) * b1*ppc*pow(3-8*sigma*pow(pc,4),2) / (1 + b1*ppc*(5 + 64*sigma*pow(pc,4)*pow(pb,4)));

        // Ra only for diagonal inlet pipes! Adjacent inlet pipes not yet implemented!
        double res_a0 = 2*res_fp + 2/(2*pi*lambda_grout) * (std::log(d_shanks/r_tube) + sigma*std::log((pow(r_BHE,4) + pow(d_shanks,4)) / (pow(r_BHE,4) - pow(d_shanks,4))));
        double res_a1 = res_a0 - 2/(2*pi*lambda_grout) * b1*ppc*pow(1+8*sigma*pow(pc,2)*pow(pb,2),2) / (1 - b1*ppc*(3 - 32*sigma*(pow(pc,2)*pow(pb,6) + pow(pc,6)*pow(pb,2))));

        // Calculate Rg and Rgg according to Bauer et al. 2011 and Pasquier and Marcotte 2012
        double res_g = n_pipes*res_b1 - res_fp; // total sheath grout thermal resistance
        double res_ar1 = (2 + pow(2,1/2))*res_g*(res_a1 - res_fp) / (res_g + res_a1 - res_fp);
        double res_ar2 = pow(2,1/2)*res_ar1;
        double res_gg1 = 2*res_g*res_ar1 / (2*res_g - res_ar1); // total core grout thermal resistance between adjacent pipes
        double res_gg2 = 2*res_g*res_ar2 / (2*res_g - res_ar2); // total core grout thermal resistance between diagonally oposed pipes

        return {res_g,res_gg1,res_gg2};
    }

    /*void solveDiffEqIntern(double aa, double bb, double Ti, double Timestep, double& Tf, double& Tbar)
    { // Trnsys function solveDiffEq with Timestep as input
        if(abs(aa) > 0.)
        { // Solution to differential is exponential
            Tf = (Ti + bb/aa)*std::exp(aa*Timestep) - bb/aa;
            Tbar = (Ti + bb/aa)/aa/Timestep*(std::exp(aa*Timestep) - 1.) - bb/aa;
        }
        else
        { //Solution to differential equation is linear
            Tf = bb*Timestep + Ti;
            Tbar = (Tf + Ti)/2.;
        }
    }*/

    double solveDiffEqIntern_ave(double aa, double bb, double Ti, double Timestep)
    { // Trnsys function solveDiffEq with Timestep as input
        if(aa != 0.)
        { // Solution to differential is exponential
            auto div = bb/aa;
            auto mul = aa*Timestep;
            //Tbar = (Ti + bb/aa)/aa/Timestep*(std::exp(aa*Timestep) - 1.) - bb/aa;
            return (Ti + div) / mul * (std::exp(mul) - 1.) - div;
        }
        //else
        // Tbar
        double Tf = bb*Timestep + Ti;
        return 0.5*(Tf + Ti);
    }

    double solveDiffEqIntern_end(double aa, double bb, double Ti, double Timestep)
    { // Trnsys function solveDiffEq with Timestep as input
        if(aa != 0.)
        { // Solution to differential is exponential
            auto div = bb/aa;
            //Tf = (Ti + bb/aa)*std::exp(aa*Timestep) - bb/aa;
            return (Ti + div)*std::exp(aa*Timestep) - div;
        }
        //else
        //Solution to differential equation is linear
        return bb*Timestep + Ti;
    }
}