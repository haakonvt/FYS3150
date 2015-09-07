#include <iostream>
#include <armadillo>
#include "DiffSolv.h"
#include "lib.h"

using namespace std;
using namespace arma;

int main()
{

    // ONE DIMENSION
    // Initial conditions
    double d = 1.0;     // Length of x-interval (and y)
    double T = 1;     // Simulation ending time
    int Nx = 10;        // Number of points in x-direction
    int Nt = 400;        // Number of time steps
    int N  = 1000;      // Total number of particles
    vec x  = zeros(N);  // Elements represents the position of the particles (for MC methods)

    int cutoff = 1000;  // Upper limit in sum (1D analytic solution)

    double dx   = d / Nx;
    double dt   = T / Nt;

    // Bin size for gaussian random walk
    double bin0 = 1.0/16;

    vec v  = zeros(Nx);
    vec us = 1-linspace(0,1,Nx);

    v(0) = 1;          // All neurotransmitters located at x = 0 when t = 0
    vec v_init = v-us; // 'Introduce' v = u - us, which give us v(endpoints) = 0
    v = v_init;        // Save for reuse (initial vector)

    // Make a solver object called 'solver'
    // Each member, a different solver, will save its solution to file.
    DiffusionSolvers solver(Nx, Nt, dx, dt, cutoff);

    // Compute solution at after some time, T:
    solver.ExplForwEuler(v);
    v = v_init; // Resets v back to initial cond.

    solver.ImplBackEuler(v);
    v = v_init; // Resets v back to initial cond.

    solver.CrankNicolson(v);

    // Analytic solution
    vec u_A = zeros(100);
    vec x_A = linspace(0,1,Nx);
    vec t_A = linspace(0,1,Nt);
    solver.AnalyticSolution(u_A, x_A, t_A);

    // Monte Carlo methods with uniform (constant) and gaussian distributed steplengths
    solver.UniStep(x);
    x = zeros(N); // Reset x to initial condition (all particles in zero)

    solver.GaussStep(x, bin0);


    // TWO DIMENSIONS
    // Initial conditions
    int Ny = Nx;
    mat u = zeros(Nx,Ny);

    // Boundary conditions u(x=0,y) = 1, all other elements = 0
    u.row(0)    = ones(Nx).t();

    solver.FTCS_2D(u);

    solver.BTCS_2D(u);





/*  ___________________________________
    TEST SECTION

    // TEST GAUSSIAN DIST
    double randR,randT;
    long seed = -1;
    vec rand = zeros(100000);
    double pi=4*atan(1.0);
    for (int k=0; k<100000;k++){
        randR = sqrt(-log(1-ran0(&seed)));
        randT = 2*pi*ran0(&seed);
        rand(k) = randR*cos(randT);
    }

    rand.save("randtest.dat",raw_ascii);

    mat test = randu(3,3);
    test.diag(0)(2) = 33;
    cout << "main diag= " << test.diag(0) << endl;
    cout << "main diag element 3= " << test.diag(0)(2) << endl;
    */

    /* CUBE test
    cube a = zeros(3,3,3);
    a.slice(2) = ones<mat>(3,3);

    cout << a.slice(2) << endl;
    a.save("simple_cube.dat",raw_ascii);*/

    /*RNG test
            long idum = -1;
    vec r1 = vec(10000);
    vec r2 = vec(10000);
    for(int i=0; i<10000;i++){
        r1(i) = randGauss(&idum);
        r2(i) = ran0(&idum);
    }
    r1.save("r1.dat",raw_ascii);
    r2.save("r2.dat",raw_ascii);*/

    return 0;
}
