#include "EulerPoissonFixedPoints.h"
#include <cmath>

EulerPoissonFixedPoints::EulerPoissonFixedPoints() {

}

EulerPoissonFixedPoints::~EulerPoissonFixedPoints() {

}

void EulerPoissonFixedPoints::find(const EulerPoisson *eqs, const double * M0, 
            double *M1, double *phi1, double *z1) {
    double g = sqrt(M0[0] * M0[0] + M0[1] * M0[1] + M0[2] * M0[2]);

    method_.init(eqs);

    double T = 1.0;
    double step = 0.001;
    int N_step = (int)ceil(T / step);

    double M1_a[3], phi1_a, z1_a;
    double M1_b[3], phi1_b, z1_b;
    double M2[3], phi2, z2;
    double M2_a[3], phi2_a, z2_a;
    double M2_b[3], phi2_b, z2_b;

    double A[2][2], invA[2][2];

    M1[0] = M0[0];
    M1[1] = M0[1];
    M1[2] = M0[2];

    double eps = 1e-12;
    double delta = 1e-6;
    double maxLen = 1e-2;

    while (true) {
        *phi1 = atan2(M1[1], M1[0]);
        *z1 = M1[2] / g;

        phi1_a = *phi1 + delta;
        z1_a = *z1;
        M1_a[0] = g * cos(phi1_a) * sqrt(1.0 - z1_a * z1_a);
        M1_a[1] = g * sin(phi1_a) * sqrt(1.0 - z1_a * z1_a);
        M1_a[2] = g * z1_a;

        phi1_b = *phi1;
        z1_b = *z1 + delta;
        M1_b[0] = g * cos(phi1_b) * sqrt(1.0 - z1_b * z1_b);
        M1_b[1] = g * sin(phi1_b) * sqrt(1.0 - z1_b * z1_b);
        M1_b[2] = g * z1_b;

        method_.map(eqs, 0.0, M1, M2, N_step, T);
        phi2 = atan2(M2[1], M2[0]);
        z2 = M2[2] / g;
        double dPhi = phi2 - *phi1;
        double dZ = z2 - *z1;
        if (sqrt(dPhi * dPhi + dZ * dZ) < eps) break;

        method_.map(eqs, 0.0, M1_a, M2_a, N_step, T);
        phi2_a = atan2(M2_a[1], M2_a[0]);
        z2_a = M2_a[2] / g;
        double dPhi_a = phi2_a - phi1_a;
        double dZ_a = z2_a - z1_a;

        method_.map(eqs, 0.0, M1_b, M2_b, N_step, T);
        phi2_b = atan2(M2_b[1], M2_b[0]);
        z2_b = M2_b[2] / g;
        double dPhi_b = phi2_b - phi1_b;
        double dZ_b = z2_b - z1_b;

        A[0][0] = (dPhi_a - dPhi) / delta;
        A[1][0] = (dZ_a - dZ) / delta;
        
        A[0][1] = (dPhi_b - dPhi) / delta;
        A[1][1] = (dZ_b - dZ) / delta;

        double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        
        invA[0][0] = A[1][1]; invA[0][1] = - A[0][1];
        invA[1][0] = - A[1][0]; invA[1][1] = A[0][0];

        double hPhi = (invA[0][0] * dPhi + invA[0][1] * dZ) / det;
        double hZ   = (invA[1][0] * dPhi + invA[1][1] * dZ) / det;
        
        double len = sqrt(hPhi * hPhi + hZ * hZ);
        if (len < maxLen) {
            *phi1 -= hPhi;
            *z1 -= hZ;
        }
        else {
            *phi1 -= hPhi / len * maxLen;
            *z1 -= hZ / len * maxLen;
        }
        M1[0] = g * cos(*phi1) * sqrt(1.0 - *z1 * *z1);
        M1[1] = g * sin(*phi1) * sqrt(1.0 - *z1 * *z1);
        M1[2] = g * *z1;
    }
}

void EulerPoissonFixedPoints::getTrace(const EulerPoisson *eqs, const double * M0, double *det, double *tr) {
    double g = sqrt(M0[0] * M0[0] + M0[1] * M0[1] + M0[2] * M0[2]);

    method_.init(eqs);

    double T = 1.0;
    double step = 0.001;
    int N_step = (int)ceil(T / step);

    double M1[3], phi1, z1;
    double M1_a[3], phi1_a, z1_a;
    double M1_b[3], phi1_b, z1_b;
    double M2[3], phi2, z2;
    double M2_a[3], phi2_a, z2_a;
    double M2_b[3], phi2_b, z2_b;

    double A[2][2];

    M1[0] = M0[0];
    M1[1] = M0[1];
    M1[2] = M0[2];

    double eps = 1e-12;
    double delta = 1e-6;

    phi1 = atan2(M1[1], M1[0]);
    z1 = M1[2] / g;

    phi1_a = phi1 + delta;
    z1_a = z1;
    M1_a[0] = g * cos(phi1_a) * sqrt(1.0 - z1_a * z1_a);
    M1_a[1] = g * sin(phi1_a) * sqrt(1.0 - z1_a * z1_a);
    M1_a[2] = g * z1_a;

    phi1_b = phi1;
    z1_b = z1 + delta;
    M1_b[0] = g * cos(phi1_b) * sqrt(1.0 - z1_b * z1_b);
    M1_b[1] = g * sin(phi1_b) * sqrt(1.0 - z1_b * z1_b);
    M1_b[2] = g * z1_b;

    method_.map(eqs, 0.0, M1, M2, N_step, T);
    phi2 = atan2(M2[1], M2[0]);
    z2 = M2[2] / g;
    double dPhi = phi2 - phi1;
    double dZ = z2 - z1;

    method_.map(eqs, 0.0, M1_a, M2_a, N_step, T);
    phi2_a = atan2(M2_a[1], M2_a[0]);
    z2_a = M2_a[2] / g;
    double dPhi_a = phi2_a - phi1_a;
    double dZ_a = z2_a - z1_a;

    method_.map(eqs, 0.0, M1_b, M2_b, N_step, T);
    phi2_b = atan2(M2_b[1], M2_b[0]);
    z2_b = M2_b[2] / g;
    double dPhi_b = phi2_b - phi1_b;
    double dZ_b = z2_b - z1_b;

    A[0][0] = (dPhi_a - dPhi) / delta;
    A[1][0] = (dZ_a - dZ) / delta;
        
    A[0][1] = (dPhi_b - dPhi) / delta;
    A[1][1] = (dZ_b - dZ) / delta;

    *det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    *tr = fabs(A[0][0] + A[1][1]);
}

