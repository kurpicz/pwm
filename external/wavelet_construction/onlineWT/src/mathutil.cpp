#include <math.h>

#include "mathutil.hpp"


double multfloor(double value, double base)
{
    if (base>0) {
        if (value>=0) {
            return floor(value/base);
        }
        else {
            return -ceil(-value/base);
        }
    }
    else if (base<0) {
        if (value<=0) {
            return ceil(value/base);
        }
        else {
            return -floor(-value/base);
        }
    }
    else { // base==0
        if (value>=0) {
            return 0;
        }
        else {
            return -INFINITY;
        }
    }
}

double multceil(double value, double base)
{
    if (base>0) {
        if (value>=0) {
            return ceil(value/base);
        }
        else {
            return -floor(-value/base);
        }
    }
    else if (base<0) {
        if (value<=0) {
            return floor(value/base);
        }
        else {
            return -ceil(-value/base);
        }
    }
    else { // base==0
        if (value<=0) {
            return 0;
        }
        else {
            return INFINITY;
        }
    }
}
