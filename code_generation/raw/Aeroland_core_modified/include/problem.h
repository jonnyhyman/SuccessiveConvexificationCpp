
#ifndef PROBLEM_H
#define PROBLEM_H

#include "ecos.h"

#ifdef __cplusplus
extern "C" // Let compliler know this is C stuff
{
#endif

#include "ecos.h"


double config[7][1];
double x_init[1][14];
double x_final[1][14];

double A_bar_parm[50][196];
double B_bar_parm[50][42];
double C_bar_parm[50][42];
double S_bar_parm[50][14];
double z_bar_parm[50][14];
double X_last_parm[50][14];
double U_last_parm[50][14];
double w_delta_sigma_parm;
double w_delta_parm;
double s_last_parm;
double w_nu_parm;


idxint q[1189];

pfloat c[11018];
pfloat h[4679];
pfloat Gpr[4880];
pfloat b[8402];
pfloat Apr[31185];

idxint Gir[4880];
idxint Gjc[11019];
idxint Air[31185];
idxint Ajc[11019];

void gather_matrices();

pwork* solver_work;

int solve();

#ifdef __cplusplus
}
#endif


#endif //PROBLEM_H

