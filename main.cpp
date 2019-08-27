#include <cstdio>
#include <cmath>
#include <valarray>
#include "Mathieu.h"
#include "RungeKutta.h"
#include <omp.h>

using namespace std;

void saveToFile(const char * fName, valarray<double> & x, valarray<double> & y, 
        valarray<valarray<double> > & T);

int main(int argc, char * argv[]) {
    double omega_min = 0.1;
    double omega_max = 2.0;
    double e_min = 0.0;
    double e_max = 2.0;

    int N_omega = 501;
    int N_e = 501;
    
    double dOmega = (omega_max - omega_min) / (N_omega - 1);
    double de = (e_max - e_min) / (N_e - 1);

    valarray<double> e(N_e);
    valarray<double> omega(N_omega);
    valarray<valarray<double> > tr(N_omega), det(N_omega);
    for (int i = 0; i < N_omega; i++) {
        tr[i].resize(N_e, 0.0);
        det[i].resize(N_omega, 0.0);
    }

    for (int i = 0; i < N_omega; i++) {
        omega[i] = omega_min + i * dOmega;
    }

    for (int j = 0; j < N_e; j++) {
        e[j] = e_min + j * de;
    }

    double T = 2.0 * M_PI;
    int N_t = 6000;
    printf("Integration Step: %lg\n", T / N_t);
    double delta = 1e-4;
    double ic1_p[2] = { delta, 0.0 }, ic1_m[2] = {- delta, 0.0 };
    double ic2_p[2] = { 0.0, delta }, ic2_m[2] = { 0.0, -delta };
   
    printf("Number of threads: ");
    int nThreads = 1;
    scanf("%d", &nThreads);
    omp_set_num_threads(nThreads);
    int count = 0;
    
    double start = omp_get_wtime();
    #pragma omp parallel
    {
        RungeKutta method;
        Mathieu eqs;
        method.init(&eqs);

        double img1_p[2], img1_m[2], img2_p[2], img2_m[2];
        double A[2][2];
        #pragma omp for
        for (int i = 0; i < N_omega; i++) {
            for (int j = 0; j < N_e; j++) {
                eqs.setParameter(0, omega[i]);
                eqs.setParameter(1, e[j]);

                method.map(&eqs, 0.0, ic1_p, img1_p, N_t, T);
                method.map(&eqs, 0.0, ic1_m, img1_m, N_t, T);
                method.map(&eqs, 0.0, ic2_p, img2_p, N_t, T);
                method.map(&eqs, 0.0, ic2_m, img2_m, N_t, T);
            
                A[0][0] = (img1_p[0] - img1_m[0]) / (2.0 * delta);
                A[1][0] = (img1_p[1] - img1_m[1]) / (2.0 * delta);
                A[0][1] = (img2_p[0] - img2_m[0]) / (2.0 * delta);
                A[1][1] = (img2_p[1] - img2_m[1]) / (2.0 * delta);
                /*A[0][0] = img1_p[0];
                A[1][0] = img1_p[1];
                A[0][1] = img2_p[0];
                A[1][1] = img2_p[1];*/

                det[i][j] = A[0][0] * A[1][1] - A[1][0] * A[0][1];
                tr[i][j] = fabs(A[0][0] + A[1][1]);
            }
            
            #pragma omp critical
            {
                count += N_e;
                printf("Ready %6d from %6d\n", count, N_omega * N_e);
            }
        }
    }
    printf("Time: %lg\n", omp_get_wtime() - start);
    saveToFile("trace.csv", omega, e, tr);
    saveToFile("det.csv", omega, e, det);

    return 0;
}


void saveToFile(const char * fName, valarray<double> & x, valarray<double> & y, 
        valarray<valarray<double> > & T) {
    FILE *f = fopen(fName, "w");
    int N_x = (int)x.size();
    int N_y = (int)y.size();

    for (int j = 0; j < N_x; j++) {
        fprintf(f, "\t%.15lg", x[j]);
    }
    fprintf(f, "\n");

    for (int i = 0; i < N_y; i++) {
        fprintf(f, "%.15lg", y[i]);
        for (int j = 0; j < N_x; j++) {
            fprintf(f, "\t%.15lg", T[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}
