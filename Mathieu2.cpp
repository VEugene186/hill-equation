#include "Mathieu2.h"
#include <cmath>

//parameters_[0] = delta
//parameters_[1] = epsilon
Mathieu2::Mathieu2() : Equation(2, 2) {

}

Mathieu2::~Mathieu2() {

}

void Mathieu2::RHS(double t, const double * q, double * dq) const {
    dq[0] = q[1];
    dq[1] = - (parameters_[0] + parameters_[1] * cos(t)) * q[0];
}
