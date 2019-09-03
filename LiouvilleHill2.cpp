#include "LiouvilleHill2.h"
#include <cmath>

//parameters_[0] = I1(0)
//parameters_[1] = I2(0)
//parameters_[2] = I3(0)
//parameters_[3] = dJ2
//parameters_[4] = Omega
//parameters_[5] = k1
LiouvilleHill2::LiouvilleHill2() : Equation(2, 6),
        I10_(parameters_[0]), I20_(parameters_[1]), I30_(parameters_[2]),
        dJ2_(parameters_[3]), Omega_(parameters_[4]), k1_(parameters_[5]) {

}

LiouvilleHill2::~LiouvilleHill2() {

}

void LiouvilleHill2::RHS(double t, const double *q, double *dq) const {
    double I1 = I10_ + pow(dJ2_ * sin(2.0 * M_PI * t), 2);
    double I2 = I20_;
    double I3 = I30_ + pow(dJ2_ * sin(2.0 * M_PI * t), 2);

    double f = 1.0 / I2 - (1.0 - k1_) / I1;
    double g = (1.0 - k1_) / I1 - 1.0 / I3;

    dq[0] = 2.0 * M_PI / Omega_ * g * q[1];
    dq[1] = 2.0 * M_PI / Omega_ * f * q[0];
}
