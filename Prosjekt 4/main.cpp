#include <iostream>
#include <armadillo>
#include "DiffSolv.h"

using namespace std;
using namespace arma;

int main()
{
    // Initial conditions
    double d = 1; // Length of x-interval
    double T = 1; // Simulation ending time
    int Nx = 10;  // Number of points in x-direction
    int Nt = 200; // Number of time steps

    double dx = d / Nx;
    double dt = T / Nt;

    vec v  = zeros(Nx);
    vec us = 1-linspace(0,1,Nx);
    v(0)   = 1;    // All neurotransmitters located at x = 0 when t = 0
    vec v_init = v-us; // 'Introduce' v = u - us, which give us v(endpoints) = 0
    v = v_init;    // Save for reuse of this initial vector

    // Compute solution at after some time, T:
    ExplForwEuler(v,Nx,Nt,dx,dt);
    v = v_init;    // Resets v back to initial cond.

    ImplBackEuler(v,Nx,Nt,dx,dt);
    v = v_init;    // Resets v back to initial cond.

    CrankNicolson(v,Nx,Nt,dx,dt);


    vec u = zeros(100);
    vec x = linspace(0,1,Nx);
    vec t = linspace(0,1,Nt);
    t = t%t; // More t-points for early values
    int cutoff = 1000;

    AnalyticSolution(u, x, t, cutoff, Nx, Nt);

    return 0;
}
