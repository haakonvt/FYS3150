#include "system.h"

using namespace std;
using namespace arma;

system::system()
{
}


mat system::RK4(mat A, vec m, double dt, int Number_of_planets){
    mat K1,K2,K3,K4;
    K1 = dt * force(A,m,      Number_of_planets);
    K2 = dt * force(A+K1/2,m, Number_of_planets);
    K3 = dt * force(A+K2/2,m, Number_of_planets);
    K4 = dt * force(A+K3,m,   Number_of_planets);
    return A + (K1 + 2*K2 + 2*K3 + K4)/6;
}

mat system::force(mat A, vec m, int Number_of_planets){
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

mat system::VER(mat A, mat pos_prev, vec m, double dt, int Number_of_planets){

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

mat system::solver(int method, double N......){
double pi = 4*atan(1.0);
double Gm = 4*pi*pi;

}
