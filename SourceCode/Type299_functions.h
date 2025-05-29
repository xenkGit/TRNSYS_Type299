#pragma once // nur 1 x kompilieren
namespace t299
{
    void message(int crntUnit, int crntType, const char errorType[], const char message[]);
    double calcNusseltNo(double re, double pr, int crntUnit, int crntType);
    std::tuple<double,double,double>calcGroutResistances(double res_fp, double lambda_grout, double lambda_ground, double r_BHE, double r_tube, double d_shanks, int n_pipes);
    //void solveDiffEqIntern(double aa, double bb, double Ti, double Timestep, double& Tf, double& Tbar);
    double solveDiffEqIntern_ave(double aa, double bb, double Ti, double Timestep);
    double solveDiffEqIntern_end(double aa, double bb, double Ti, double Timestep);
}