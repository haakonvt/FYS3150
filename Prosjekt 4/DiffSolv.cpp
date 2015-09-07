#include "DiffSolv.h"


vec TriDiag(vec aa, vec ab, vec ac, vec b, int Nx){

    vec v = zeros(Nx); // Solution vector

    //aa(0) = 0;
    //ac(Nx-1) = 0;

    // Forward substitution
    ac(0) = ac(0)/ab(0);
    b(0)  = b(0)/ab(0);

    for (int i=1; i<=Nx-2; i++){
        ac(i) = ac(i) / (ab(i) - aa(i)*ac(i-1));
        b(i)  = (b(i) - aa(i)*b(i-1)) / (ab(i) - aa(i)*ac(i-1));
    }

    b(Nx-1) = (b(Nx-1)-aa(Nx-1)*b(Nx-2)) / (ab(Nx-1)-aa(Nx-1)*ac(Nx-2));

    // Backwards substitution gives the solutions
    v(Nx-1) = b(Nx-1);

    for (int i=Nx-1; i>=1; i--){
        v(i-1) = b(i-1) - ac(i-1)*v(i);
    }
    return v;

}


void ExplForwEuler(vec &v, int Nx, int Nt, double dx, double dt){
    vec v_next = zeros(Nx);
    mat sol    = zeros(Nx,Nt);
    vec us     = 1-linspace(0,1,Nx);

    sol.col(0) = v+us; // Save real solution u = v + us

    double dtxx = dt/(dx*dx);

    for(int j=1; j<Nt; j++){
        for(int i=1; i<Nx-1; i++){
            v_next(i) = dtxx*v(i-1) + (1-2*dtxx)*v(i) + dtxx*v(i+1);
        }
        v = v_next;

        sol.col(j) = v+us;
    }
    sol.save("ExplForwEuler.dat",raw_ascii);
}


void ImplBackEuler(vec &v, int Nx, int Nt, double dx, double dt){
    double dtxx = dt/(dx*dx);
    vec a = -dtxx*ones(Nx);
    vec b = (1 + 2*dtxx)*ones(Nx);
    vec c = a;

    mat sol    = zeros(Nx,Nt);
    vec us     = 1-linspace(0,1,Nx);
    sol.col(0) = v+us; // Save real solution u = v + us

    for(int j=1; j<Nt; j++){
        v = TriDiag(a,b,c,v,Nx);

        sol.col(j) = v+us;
    }

    sol.save("ImplBackEuler.dat",raw_ascii);
}


void CrankNicolson(vec &v, int Nx, int Nt, double dx, double dt){
    double d1 = dt/(dx*dx);
    double d2 = 2-2*d1;
    double d3 = 2+2*d1;

    vec a = -d1*ones(Nx);
    vec b = d3*ones(Nx);
    vec c = a;

    vec vv     = v;
    mat sol    = zeros(Nx,Nt);
    vec us     = 1-linspace(0,1,Nx);
    sol.col(0) = v+us; // Save real solution u = v + us

    for (int j=1; j<Nt; j++){
        for (int i=1; i<Nx-1; i++){
            vv(i) = d1*v(i-1) + d2*v(i) + d1*v(i+1);
            cout << "j=" << j << ", i=" << i << endl;

        }
        v = TriDiag(a, b, c, vv, Nx);
        sol.col(j) = v+us;
    }
    sol.save("CrankNicolson.dat",raw_ascii);

}


void AnalyticSolution(vec u, vec x, vec t, int cutoff, int Nx, int Nt){
    mat sol = mat(Nx,Nt);
    double pi = 4*atan(1.0);

    for(int i=0; i<Nt; i++){
        u = 1-x;
        for(int n=1; n<cutoff; n++){
            u += -2*sin(n*pi*x)*exp(-n*n*pi*pi*t(i))/(pi*n);
        }
        sol.col(i) = u;
    }

    sol.save("u_analytic.dat", raw_ascii);

}
