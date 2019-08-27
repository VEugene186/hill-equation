#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "Equation.h"

class RungeKutta {
public:
    RungeKutta();
    ~RungeKutta();
    void init(const Equation *eq);
    void makeStep(const Equation *eq, double t0, const double *q0, double *q1, double dt);
    void map(const Equation *eq, double t0, const double *q0, double *q1, int N_steps, double T);
private:
    int dim_;
    double *k1_, *k2_, *k3_, *k4_, *q_, *qTmp_;


};

#endif
