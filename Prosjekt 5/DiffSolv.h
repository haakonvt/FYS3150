#ifndef DIFFSOLV_H
#define DIFFSOLV_H

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class DiffusionSolvers{
    int Nx;
    int Nt;
    double dx;
    double dt;
    int cutoff;
    double pi;
    double dtxx;
    double dtxx2;
    int t,i,j,k,n; // Indices for a bunch of loops

public:
    DiffusionSolvers(int, int, double, double, int);
    void FTCS_2D(mat);
    void BTCS_2D(mat);
    vec TriDiag(vec, vec, vec, vec, int);
    void ExplForwEuler(vec &);
    void ImplBackEuler(vec &);
    void CrankNicolson(vec &);
    void AnalyticSolution(vec, vec, vec);
    void UniStep(vec);
    void GaussStep(vec, double bin0);

};

#endif // DIFFSOLV_H
