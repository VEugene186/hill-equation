#ifndef LIOUVILLEHILL_H
#define LIOUVILLEHILL_H

#include "Equation.h"

class LiouvilleHill : public Equation {
public:
    LiouvilleHill();
    ~LiouvilleHill();

    virtual void RHS(double t, const double *q, double *dq) const;
private:
    double &I10_, &I20_, &I30_, &dJ2_, &Omega_;
};

#endif
