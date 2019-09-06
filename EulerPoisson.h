#ifndef EULERPOISSON_H
#define EULERPOISSON_H

#include "Equation.h"

class EulerPoisson : public Equation {
public:
    EulerPoisson();
    ~EulerPoisson();

    virtual void RHS(double t, const double *q, double *dq) const;

private:
    double &I10_, &I20_, &I30_, &k1_, &k2_, &k3_, &Omega_, &dJ2_;

};

#endif
