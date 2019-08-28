#include "LiouvilleHill.h"
#include <cmath>

//parameters_[0] = I1(0)
//parameters_[1] = I2(0)
//parameters_[2] = I3(0)
//parameters_[3] = dJ2
//parameters_[4] = Omega
LiouvilleHill::LiouvilleHill() : Equation(2, 5),
        I10_(parameters_[0]), I20_(parameters_[1]), I30_(parameters_[2]),
        dJ2_(parameters_[3]), Omega_(parameters_[4]) {

}

LiouvilleHill::~LiouvilleHill() {

}

void LiouvilleHill::RHS(double t, const double *q, double *dq) const {
    double I1 = I10_ + pow(dJ2_ * sin(2.0 * M_PI * t), 2);
    double I2 = I20_;
    double I3 = I30_ + pow(dJ2_ * sin(2.0 * M_PI * t), 2);

    double dI1 = Omega_ * dJ2_ * dJ2_ * sin(4.0 * M_PI * t);
    double d2I1 = 2.0 * Omega_ * Omega_ * dJ2_ * dJ2_ * cos(4.0 * M_PI * t);
    double dI2 = 0.0;
    double d2I2 = 0.0;
    double dI3 = dI1;
    double d2I3 = d2I1;

    double f = 1.0 / I2 - 1.0 / I1;
    double g = 1.0 / I1 - 1.0 / I3;

    double dg = dI3 / (I3 * I3) - dI1 / (I1 * I1);
    double d2g = 2.0 * (dI1 * dI1 / (I1 * I1 * I1) - dI3 * dI3 / (I3 * I3 * I3)) + 
        d2I3 / (I3 * I3) - d2I1 / (I1 * I1);

    double p = (2.0 * d2g * g - 3.0 * dg * dg) / (4.0 * g * g) - f * g;

    dq[0] = q[1];
    dq[1] = -p * pow(2.0 * M_PI / Omega_, 2) * q[0];
}
