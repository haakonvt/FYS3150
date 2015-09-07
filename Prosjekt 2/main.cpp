#include <iostream>
#include <armadillo>
#include <eig_jac.cpp>
#include <init_one_elec.cpp>
#include <init_two_elec.cpp>

using namespace std;
using namespace arma;

int main()
{
    int N_step     = 200;
    int counter    = 0;
    double rho_max = 5;
    double wr      = 1;

    mat A = zeros(N_step-1,N_step-1);
    vec arma_eigval;
    mat eigvec;
    mat eigv1 = zeros(N_step-1,2);

    // Initial conditions and matrix configuration (one or two electrons in HO-potential) **COMMENT/UNCOMMENT**
    //init_one(N_step,A,rho_max,eigv1);
    init_two(N_step,A,rho_max,eigv1,wr);

    // Use armadillos own algorithm for getting the eig.values and eig.vectors
    // (A must be symmetric!)
    clock_t start1, finish1; // Timing computations
    start1 = clock();

    eig_sym(arma_eigval, eigvec, A);

    finish1 = clock();
    double comp_time1 = ((finish1-start1)/(double) CLOCKS_PER_SEC);

    // Use my own general solver, Jacobi rotations

    double tol = 1E-8; // Tolerance of how small off-diag. elements must be for algo. to stop
    mat C = abs(A);    // Take absolute value of all elements (used to find the max value)
    C.diag(0) = zeros(N_step-1); // Sets diag. elements equal to zero (to remove from max value search)

    // Search through elements to find max value and changes row & col index accordingly
    uword  k; // Index of row maximum
    uword  l; // Index of col maximum
    double max_val = C.max(k,l);

    //JACOBIS ROTATIONS
    clock_t start2, finish2;
    start2 = clock();

    jacobi_method(max_val,A,C,counter,tol,l,k,N_step);

    // Eigenvalues must be sorted from smallest to largest
    vec jacobi_eig_val = sort(A.diag(0));

    finish2 = clock();
    double comp_time2 = ((finish2-start2)/(double) CLOCKS_PER_SEC);

    cout << "Eig1 (jacobi   )= " << jacobi_eig_val(0) << ". Eig2= " << jacobi_eig_val(1) << ". Eig3= " << jacobi_eig_val(2) << endl;
    cout << "Eig1 (armadillo)= " << arma_eigval(0) << ". Eig2= " << arma_eigval(1) << ". Eig3= " << arma_eigval(2) << endl;
    //cout << "Percent error= " << 100*(3-jacobi_eig_val(0))/3 << ".P.err= " << 100*(7-jacobi_eig_val(1))/7 << ". P.err= " << 100*(11-jacobi_eig_val(2))/11 << endl;
    cout << "Number of runs through the program loop: " << counter << endl;
    cout << "Time taken (sec): Arma: " << comp_time1 << " , Jacobi: " << comp_time2 << endl;//

    // Save the three first eigenvectors to files
    eigv1.col(1) = eigvec.col(0);
    eigv1.save("eig1.dat", raw_ascii);
    eigv1.col(1) = eigvec.col(1);
    eigv1.save("eig2.dat", raw_ascii);
    eigv1.col(1) = eigvec.col(2);
    eigv1.save("eig3.dat", raw_ascii);

    return 0;
}


