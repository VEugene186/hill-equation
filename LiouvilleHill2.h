#ifndef LIOUVILLEHILL2_H
#define LIOUVILLEHILL2_H

#include "Equation.h"

class LiouvilleHill2 : public Equation {
public:
    LiouvilleHill2();
    ~LiouvilleHill2();

    virtual void RHS(double t, const double *q, double *dq) const;
private:
    double &I10_, &I20_, &I30_, &dJ2_, &Omega_, &k1_;
};

#endif
