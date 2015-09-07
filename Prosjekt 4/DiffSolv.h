#ifndef DIFFSOLV_H
#define DIFFSOLV_H

#endif // DIFFSOLV_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

vec TriDiag(vec aa, vec ab, vec ac, vec b, int Nx);
void ExplForwEuler(vec &v, int Nx, int Nt, double dx, double dt);
void ImplBackEuler(vec &v, int Nx, int Nt, double dx, double dt);
void CrankNicolson(vec &v, int Nx, int Nt, double dx, double dt);
void AnalyticSolution(vec u, vec x, vec t, int cutoff, int Nx, int Nt);
