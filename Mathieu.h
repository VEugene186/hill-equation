#ifndef MATHIEU_H
#define MATHIEU_H

#include "Equation.h"

class Mathieu : public Equation {
public:
    Mathieu();
    ~Mathieu();

   virtual void RHS(double t, const double * q, double * dq) const;
};

#endif
