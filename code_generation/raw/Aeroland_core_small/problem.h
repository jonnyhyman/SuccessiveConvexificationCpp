
#ifndef PROBLEM_H
#define PROBLEM_H

#include "ecos.h"


double x_init[1][14];
double x_final[1][14];
double config[7][1];
double A_bar_parm[5][196];
double B_bar_parm[5][42];
double C_bar_parm[5][42];
double S_bar_parm[5][14];
double z_bar_parm[5][14];
double X_last_parm[5][14];
double U_last_parm[5][14];
double s_last_parm[1][1];
double w_delta_parm[1][1];
double w_nu_parm[1][1];
double w_delta_sigma_parm[1][1];



idxint q[109];


pfloat c[1073];
pfloat h[449];
pfloat Gpr[470];
pfloat b[842];
pfloat Apr[2835];


idxint Gir[470];
idxint Gjc[1074];
idxint Air[2835];
idxint Ajc[1074];



void gather_matrices();

pwork* solver_work;

int solve();

void clean();

#endif //PROBLEM_H