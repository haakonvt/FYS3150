#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;


void init_one(int N_step, mat &A, double rho_max, mat &eigv1){
    double rho_min  = 0;
    double h        = (rho_max-rho_min)/N_step;


    // Set up vector rho and matrix A
    vec rho = rho_min + linspace(0,N_step,N_step+1)*h;  // rho = rho_min + i* h, i=0,1,2...N_step
    vec V   = rho%rho;                                  // Elementwise multiplication

    eigv1.col(0) = rho.subvec(1,N_step-1); // Saves for later output to file

    double ihh =  1/(h*h);
    A.diag(0)  =  2*ihh + V.subvec(1,N_step-1);
    A.diag(1)  = -ihh*ones(N_step-2);
    A.diag(-1) = -ihh*ones(N_step-2);

    return;
}
