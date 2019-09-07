#ifndef EULERPOISSONFIXEDPOINS_H
#define EULERPOISSONFIXEDPOINS_H

#include "EulerPoisson.h"
#include "RungeKutta.h"

class EulerPoissonFixedPoints {
public:
    EulerPoissonFixedPoints();
    ~EulerPoissonFixedPoints();

    void find(const EulerPoisson *eqs, const double * M0, double *M1, double *phi1, double *z1);
    void getTrace(const EulerPoisson *eqs, const double * M0, double *det, double *tr);
private:
    RungeKutta method_;
};

#endif
