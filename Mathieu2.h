#ifndef MATHIEU2_H
#define MATHIEU2_H

#include "Equation.h"

//This class duplicate class Mathieu.
//There are differences in definition of the parameters
//I try compare my calculations with results of a paper
//I.Kovacic, R.Rand, S.M.Sah Mathieu's Equation and Iis Generalizations:...
class Mathieu2 : public Equation {
public:
    Mathieu2();
    ~Mathieu2();

   virtual void RHS(double t, const double * q, double * dq) const;
};

#endif
