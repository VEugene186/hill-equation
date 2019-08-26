#include "Mathieu.h"
#include <cmath>

//parameters_[0] = omega
//parameters_[1] = epsilon
Mathieu::Mathieu() : Equation(2, 2) {

}

Mathieu::~Mathieu() {

}

void Mathieu::RHS(double t, const double * q, double * dq) const {
    dq[0] = q[1];
    dq[1] = - pow(parameters_[0], 2) * (1 + parameters_[1] * cos(t)) * q[0];
}
