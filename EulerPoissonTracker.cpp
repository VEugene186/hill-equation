#include "EulerPoissonTracker.h"
#include <cmath>
#include <cstdio>

EulerPoissonTracker::EulerPoissonTracker() {

}

EulerPoissonTracker::~EulerPoissonTracker() {

}

void EulerPoissonTracker::track(EulerPoisson * eqs, const double *M0, 
        double dJ2_min, double dJ2_max, int N_J, valarray<double*> &points, valarray<double> & dJ2) {

    double M1[3], phi1, z1;
    double g = sqrt(M0[0] * M0[0] + M0[1] * M0[1] + M0[2] * M0[2]);

    double step = (dJ2_max - dJ2_min) / (N_J - 1);
    points.resize(N_J);
    dJ2.resize(N_J);
    points[0] = new double[5];
    
    dJ2[0] = dJ2_min;
    eqs->setParameter(7, dJ2_min);
    epfp_.find(eqs, M0, points[0], points[0] + 3, points[0] + 4);
    /*int i = 0;
    printf("%.15lg\t%.15lg\t%.15lg\t%.15lg\t%.15lg\n", 
            points[i][0], points[i][1], points[i][2], 
            points[i][3], points[i][4]);*/
    
    double icM[5];
    for (int i = 1; i < N_J; i++) {
        points[i] = new double[5];
        dJ2[i] = dJ2_min + i * step;
        eqs->setParameter(7, dJ2[i]);
        if (i > 1) {
            icM[3] = 2.0 * points[i - 1][3] - points[i - 2][3];
            icM[4] = 2.0 * points[i - 1][4] - points[i - 2][4];
            icM[0] = g * cos(icM[3]) * sqrt(1.0 - icM[4] * icM[4]);
            icM[1] = g * sin(icM[3]) * sqrt(1.0 - icM[4] * icM[4]);
            icM[2] = g * icM[4];
        }
        else {
            icM[0] = points[0][0];
            icM[1] = points[0][1];
            icM[2] = points[0][2];
            icM[3] = points[0][3];
            icM[4] = points[0][4];
        }

        epfp_.find(eqs, icM, points[i], points[i] + 3, points[i] + 4);
        /*printf("%8lg | ", dJ2_min + i * step);
        printf("%.15lg\t%.15lg\t%.15lg\t%.15lg\t%.15lg\n", 
                points[i][0], points[i][1], points[i][2], 
                points[i][3], points[i][4]);*/
    }
}
