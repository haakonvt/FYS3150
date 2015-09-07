#include <iostream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

void jacobi_method (double max_val,mat &A, mat &C, int &counter, double tol,uword l,uword k,int N_step){

    double t,c,s,tau;
    double akk, all, aik, ail;

    while(max_val>tol){
        tau = (A(l,l)-A(k,k))/(2*A(k,l));

        // Compute t = tan(theta) two equal ways to minimize errors dep. if tau is neg or pos
        if(tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else {
            t = -1.0/( -tau + sqrt(1.0 + tau*tau));
        }

        c = 1.0/sqrt(1+ t*t);
        s = c*t;

        // Update matrix A with new computed values (and zeros!)

        akk = A(k,k);
        all = A(l,l);

        A(k,k) = c*c*akk - 2*c*s*A(k,l) + s*s*all;
        A(l,l) = s*s*akk + 2*c*s*A(k,l) + c*c*all;
        A(k,l) = 0; // Make sure values are exactly zero
        A(l,k) = 0;

        // Update the rest of the affected elements
        for(int i=0; i<N_step-1; i++){
            if(i != k && i != l ) {  // Must skip already updated elements
                aik = A(i,k);
                ail = A(i,l);
                A(i,k) = c*aik - s*ail;
                A(k,i) = A(i,k);
                A(i,l) = c*ail + s*aik;
                A(l,i) = A(i,l);
            }
        }

        // New search for the _now_ largest element in A
        C = abs(A);
        C.diag(0) = zeros(N_step-1);
        max_val = C.max(k,l);

        // Count number of runs through the loop
        counter++;

    }
}
