#ifndef SYSTEM_H
#define SYSTEM_H
#include <armadillo>

class system
{
public:
    system();

    arma::mat force();

    arma::mat RK4(arma::mat A, arma::vec m, double dt, int Number_of_planets);

    arma::mat VER(mat A, mat pos_prev, vec m, double dt, int Number_of_planets);
};

#endif // SYSTEM_H
