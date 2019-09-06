#ifndef EULERPOISSONTRACKER_H
#define EULERPOISSONTRACKER_H

#include "EulerPoissonFixedPoints.h"
#include <valarray>
using namespace std;

class EulerPoissonTracker {
public:
    EulerPoissonTracker();
    ~EulerPoissonTracker();

    void track(EulerPoisson * eqs, const double *M0, double dJ2_min, double dJ2_max,
                int N_J, valarray<double*> &points, valarray<double> & dJ2);
private:
    EulerPoissonFixedPoints epfp_;
};

#endif
