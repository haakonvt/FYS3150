#ifndef PLANET_H
#define PLANET_H
#include <armadillo>

class planet
{
public:
    double m,x,vy;
    planet(double x, double vy, double m);
};

#endif // PLANET_H
