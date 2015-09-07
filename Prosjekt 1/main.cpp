#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;

int main()
{
    // Measure time taken solving the problem
    clock_t start, finish;
    start = clock();

    int n = 1000;
    double h = 1.0/(n+1);
    vec x = linspace(h,1-h,n);

    // Produce b-vector with elements: b_i = h^2 * f(x_i)
    double h2 = h*h;
    vec b = zeros(n);
    vec f = 100*exp(-10*x);
    b.fill(h2);
    b = b%f; // Element-wise multiplication
    vec d = b;   // Saves original RHS for other solution methods

    // Write tridiagonal matrix as vectors
    vec aa = -1*ones(n);
    vec ac = aa;
    vec ab = 2*ones(n);

    // Dirichlet bcs
    aa(0) = 0;
    ac(n-1) = 0;

    // Forward substitution
    // See my report for the derivation
    ac(0) = ac(0)/ab(0);
    b(0)  = b(0)/ab(0);

    for (int i=1; i<=n-2; i++){
        ac(i)   = ac(i) / (ab(i) - aa(i)*ac(i-1));
        b(i)  = (b(i) - aa(i)*b(i-1)) / (ab(i) - aa(i)*ac(i-1));
    }

    b(n-1) = (b(n-1)-aa(n-1)*b(n-2)) / (ab(n-1)-aa(n-1)*ac(n-2));

    // Backwards substitution gives the solutions
    vec v = zeros(n); // Solution vector

    v(n-1) = b(n-1);

    for (int i=n-1; i>=1; i--){
        v(i-1) = b(i-1) - ac(i-1)*v(i);
    }

    // Stop timing of computation
    finish = clock();
    double comp_time = ((finish-start)/(double) CLOCKS_PER_SEC);

   // Write solution and x-values to file
   // ..for plotting in external program
    mat to_file = zeros(n+2,2);
    for (int j=1; j<n+1; j++){
       to_file(j,0) = x(j-1);
       to_file(j,1) = v(j-1);
    }
    to_file(n+1,0) = 1;
    to_file.save("sol_vec_n20.dat", raw_ascii);

    // Compute the relative error between numerical and exact solution.
    vec e = zeros(n);
    vec u = 1-(1-exp(-10))*x-exp(-10*x);
    e = log10(abs((v-u)/u));
    e.save("rel_error_n20.dat",raw_ascii);
    cout << "Maximum relative error for n= " << n << " is " << max(e) << endl;

    // Compare time taken to solve problem with standard Gaussian elimination and LU decomposition.
    // Setup of tridiagonal matrix (now actually as a matrix)
    mat A = 2*eye(n,n);             // Identity matrix times 2 fills diagonal "i=j"
    A.diag(-1) = -1*ones(n-1);      // Fills diagonals above and beneath main diag. with "-1"s
    A.diag(+1) = -1*ones(n-1);

    clock_t start2, finish2;
    start2 = clock();

    // Solve Ax=b with inbuilt Gaussian elimination
    vec vGauss = solve(A,d);

    // Stop clock
    finish2 = clock();
    double comp_time_gauss = ((finish2-start2)/(double) CLOCKS_PER_SEC);

    // Last one out, LU decomposition
    // P stores possible row interchanges, so that A = trans(P)*L*U
    mat P,L,U;
    lu(L, U, P, A);

    // Assume LU decomp. is known
    clock_t start3, finish3;

    vec y   = solve(L,d);
    vec vLU = solve(U,y);

    finish3 = clock();
    double comp_time_LU = ((finish3-start3)/(double) CLOCKS_PER_SEC);

    /*
    Check answers for our three algorithms
    cout << v << endl;
    cout << vGauss << endl;
    cout << vLU << endl;
    */

    cout << "For n= " << n << ". Time taken:" << endl;
    cout << "Own solver: " << comp_time << " seconds." << endl;
    cout << "Armadillo LU decomp. solve: " << comp_time_LU << " seconds." << endl;
    cout << "Armadillo solve(A,b): " << comp_time_gauss << " seconds." << endl;


    return 0;
}

