#include <iostream>
#include <armadillo>
#include <planet.h>
#include <system.h>

using namespace std;
using namespace arma;

int main()
{
    int Number_of_planets = 10; // (Pluto & the sun included)
    double year  = 500;
    int steps_yr = 200 ;        // Steps per year
    int N        = year*steps_yr;
    double dt    = year/N;

    // Columns will represent each timestep
    // Rows    will represent each planet ("object")
    mat A  = zeros(Number_of_planets,4);
    mat x  = zeros(Number_of_planets,N);
    mat y  = zeros(Number_of_planets,N);
    vec m  = zeros(Number_of_planets);
    vec vy = zeros(Number_of_planets);


    // Average values used for distance and velocity (length = [au], time = [year])
    // (will give circular orbits)    // planet some_planet(x,vy,mass)
    planet Sun    (0,        0.00206169, 19.89);
    planet Mercury(-0.387, -10.098,      0.00000330);
    planet Venus  (0.722,    7.378,      0.0000487);
    planet Earth  (-1,      -6.282,      0.00005974);
    planet Mars   (1.52,     5.0963,     0.00000642);
    planet Jupiter(-5.2,    -2.76,       0.01899);
    planet Saturn (9.58,     2.023,      0.00568);
    planet Uranus (-19.2,   -1.433,      0.000866);
    planet Neptune(30.1,     1.138,      0.00103);
    planet Pluto  (39.5,     0.991,      0.0000001309);

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
    xx.col(0) = x.col(0);               // Positions of all planets/objects
    A = zeros(Number_of_planets,4);   // Clear A for new use
    A.col(0) = xx.col(0);               // Setup of A for first step
    A.col(3) = vy;

    // First step using RK4
    A = RK4(A, m, dt, Number_of_planets);
    xx.col(1) = A.col(0);
    yy.col(1) = A.col(1);

    // Matrix for position at previous timestep
    mat pos_prev    = zeros(Number_of_planets,4);
    pos_prev.col(0) = xx.col(0);
    pos_prev.col(1) = yy.col(0);

    /* Matrix  containing energy and ang. momentum
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
    x.save("x.dat",  raw_ascii);
    y.save("y.dat",  raw_ascii);
    xx.save("xx.dat",raw_ascii);
    yy.save("yy.dat",raw_ascii);
    //E.save("energy.dat",raw_ascii);
    return 0;
}















