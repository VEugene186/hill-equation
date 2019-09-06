#include "EulerPoisson.h"
#include <cmath>

EulerPoisson::EulerPoisson() : Equation(3, 8),
        I10_(parameters_[0]), I20_(parameters_[1]), I30_(parameters_[2]), 
        k1_(parameters_[3]), k2_(parameters_[4]), k3_(parameters_[5]),
        Omega_(parameters_[6]), dJ2_(parameters_[7])
{

}

EulerPoisson::~EulerPoisson() {

}

void EulerPoisson::RHS(double t, const double *q, double *dq) const {
    double I1 = I10_ + pow(dJ2_ * sin(2.0 * M_PI * t), 2);
    double I2 = I20_;
    double I3 = I30_ + pow(dJ2_ * sin(2.0 * M_PI * t), 2);

    double M1 = q[0];
    double M2 = q[1];
    double M3 = q[2];
    double omega1 = (q[0] - k1_) / I1;
    double omega2 = (q[1] - k2_) / I2;
    double omega3 = (q[2] - k3_) / I3;

    double k = 2.0 * M_PI / Omega_;

    dq[0] = k * (M2 * omega3 - M3 * omega2);
    dq[1] = k * (M3 * omega1 - M1 * omega3);
    dq[2] = k * (M1 * omega2 - M2 * omega1);
}

