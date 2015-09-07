#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

double pi = 4*atan(1.0);
double Gm = 4*pi*pi;

mat force(mat A, vec m, int Number_of_planets){
    mat dAdt = zeros(Number_of_planets,4);
    dAdt.col(0) = A.col(2);
    dAdt.col(1) = A.col(3);
    vec x = A.col(0);
    vec y = A.col(1);
    for (int i=0; i<Number_of_planets; i++){
        for (int j=0; j<Number_of_planets; j++){
            if (i != j){
                double posx = x(i) - x(j);
                double posy = y(i) - y(j);
                double posr = sqrt(posx*posx + posy*posy);
                // Sum up the contribution of forces from all other bodies:
                dAdt(i,2) += -Gm*m(j)*posx/(m(0)*posr*posr*posr);
                dAdt(i,3) += -Gm*m(j)*posy/(m(0)*posr*posr*posr);
            }
        }
    }
    return dAdt;
}
mat RK4(mat A, vec m, double dt, int Number_of_planets){
    mat K1,K2,K3,K4;
    K1 = dt * force(A,m, Number_of_planets);
    K2 = dt * force(A+K1/2,m, Number_of_planets);
    K3 = dt * force(A+K2/2,m, Number_of_planets);
    K4 = dt * force(A+K3,m, Number_of_planets);
    return A + (K1 + 2*K2 + 2*K3 + K4)/6;
}
mat VER(mat A, mat pos_prev, vec m, double dt, int Number_of_planets){
    vec x = A.col(0);
    vec y = A.col(1);
    vec accelx = zeros(Number_of_planets);
    vec accely = zeros(Number_of_planets);
    for (int i=0; i<Number_of_planets; i++){
        for (int j=0; j<Number_of_planets; j++){
            if (i != j){
                double posx = x(i) - x(j);
                double posy = y(i) - y(j);
                double posr = sqrt(posx*posx + posy*posy);
                // Sum up the contribution of forces from all other bodies:
                accelx(i) += -Gm*m(j)*posx/(m(0)*posr*posr*posr);
                accely(i) += -Gm*m(j)*posy/(m(0)*posr*posr*posr);
            }
        }
    }
    A.col(0) = 2*x - pos_prev.col(0) + dt*dt*accelx;
    A.col(1) = 2*y - pos_prev.col(1) + dt*dt*accely;
    return A;
}
int main()
{
    int Number_of_planets = 10; // (Pluto & the sun included)
    double year = 500;
    int steps_yr = 200 ; // Steps per year
    int N = year*steps_yr;
    double dt = year/N;
    // Columns will represent each timestep
    // Rows will represent each planet ("object")
    mat A = zeros(Number_of_planets,4);
    mat x = zeros(Number_of_planets,N);
    mat y = zeros(Number_of_planets,N);
    vec m = zeros(Number_of_planets);
    vec vy = zeros(Number_of_planets);

    // Average values used for distance and velocity (length = [au], time = [year])
    // (will give circular orbits)
    /* Three body problem
x(0,0) = 0; vy(0) = -2.6348e-05; m(0) = 19.89; // Sun
x(1,0) = 1.0; vy(1) = sqrt(1.95)*6.282; m(1) = 0.00005974; // Earth
x(2,0) = -5.20; vy(2) = -2.76; m(2) = 1000*0.01899; // Jupiter*/

    x(0,0) = 0; vy(0) = 0; m(0) = 19.89; // Sun
    x(1,0) = -0.387; vy(1) = -10.098; m(1) = 0.00000330; // Mercury
    x(2,0) = 0.722; vy(2) = 7.378; m(2) = 0.0000487; // Venus
    x(3,0) = -1.0; vy(3) = -6.282; m(3) = 0.00005974; // Earth
    x(4,0) = 1.52; vy(4) = 5.0963; m(4) = 0.00000642; // Mars
    x(5,0) = -5.20; vy(5) = -2.76; m(5) = 0.01899; // Jupiter
    x(6,0) = 9.58; vy(6) = 2.023; m(6) = 0.00568; // Saturn
    x(7,0) = -19.2; vy(7) = -1.433; m(7) = 0.000866; // Uranus
    x(8,0) = 30.1; vy(8) = 1.138; m(8) = 0.00103; // Neptune
    x(9,0) = 39.5; vy(9) = 0.991; m(9) = 0.0000001309; // Pluto

    // Net momentum must be zero to avoid movement of entire system
    vy(0) = -sum(m%vy)/m(0);

    // Setup copies for Verlet
    mat xx = zeros(Number_of_planets,N);
    mat yy = zeros(Number_of_planets,N);
    vec axx = zeros(Number_of_planets);
    vec ayy = zeros(Number_of_planets);

    // Initial conditions for RungeKutta4
    A.col(0) = x.col(0);
    A.col(3) = vy;

    // RK4 loop
    for (int i=0; i<N; i++){
        x.col(i) = A.col(0);
        y.col(i) = A.col(1);
        // Sun stays at x=y=0
        //A.row(0) = zeros(1,4);
        A = RK4(A, m, dt, Number_of_planets);
    }

    // Initial conditions for Verlet
    xx.col(0) = x.col(0); // Positions of all planets/objects
    A = zeros(Number_of_planets,4); // Clear A for new use
    A.col(0) = xx.col(0); // Setup of A for first step
    A.col(3) = vy;

    // First step using RK4
    A = RK4(A, m, dt, Number_of_planets);
    xx.col(1) = A.col(0);
    yy.col(1) = A.col(1);

    // Matrix for position at previous timestep
    mat pos_prev = zeros(Number_of_planets,4);
    pos_prev.col(0) = xx.col(0);
    pos_prev.col(1) = yy.col(0);

    /* Matrix containing energy and ang. momentum
mat E = zeros(3,N-1);
double posr,posx,posy,velx,vely,velv;*/

    for (int i=1; i<N-1; i++){
        A = VER(A, pos_prev, m, dt, Number_of_planets);
        /* Energy & ang. mom.
posx = A(1,0);
posy = A(1,1);
velx = (posx-pos_prev(1,0))/dt;
vely = (posy-pos_prev(1,1))/dt;
velv = sqrt(velx*velx + vely*vely);
posr = sqrt(posx*posx + posy*posy);
E(0,i-1) = -Gm*m(1)/posr;
E(1,i-1) = 0.5*m(1)*velv*velv;
E(2,i-1) = posr*m(1)*velv;*/
        xx.col(i+1) = A.col(0);
        yy.col(i+1) = A.col(1);
        pos_prev.col(0) = xx.col(i);
        pos_prev.col(1) = yy.col(i);
    }

    // Save data from simulation
    x.save("x.dat", raw_ascii);
    y.save("y.dat", raw_ascii);
    xx.save("xx.dat",raw_ascii);
    yy.save("yy.dat",raw_ascii);
    //E.save("energy.dat",raw_ascii);
    return 0;
}
