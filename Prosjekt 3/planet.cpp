#include "planet.h"

planet::planet(double x, double vy, double mass)
{
    posx = x;
    posy = 0;
    velx = 0;
    vely = vy;
    m    = mass;

}
