#include <cstdio>
#include <cmath>
#include <valarray>
#include <vector>
#include <string>
#include "Mathieu.h"
#include "Mathieu2.h"
#include "LiouvilleHill.h"
#include "LiouvilleHill2.h"
#include "RungeKutta.h"
#include <omp.h>

using namespace std;

void saveToFile(const char * fName, valarray<double> & x, valarray<double> & y, 
        valarray<valarray<double> > & T);

void startMenuLoop();

void multitreadingOptions();
void clearStream();

void case1();
void case2();
void case3();
void case4();

int main(int argc, char * argv[]) {
    startMenuLoop();
    
    return 0;

    printf("Number of threads: ");
    int nThreads = 1;
    scanf("%d", &nThreads);
    omp_set_num_threads(nThreads);

    printf("1 - Mathieu (Arnold)\n"
           "2 - Mathieu (Kovacic)\n"
           "3 - LiouvillHill\n"
           "4 - LiouvillHill2\n");
    int selection = -1;
    scanf("%d", &selection);
    switch (selection) {
    case 1:
        case1();
        break;
    case 2:
        case2();
        break;
    case 3:
        case3();
        break;
    case 4:
        case4();
        break;
    default:
        printf("Unknown selection\n");
    }


    return 0;
}

void clearStream() {
    while (getchar() != '\n') {}
}

void multitreadingOptions() {
    int nThreads = omp_get_max_threads();
    printf("Number of thread [%d]: ", nThreads);
    scanf("%d", &nThreads);
    omp_set_num_threads(nThreads);
}

void startMenuLoop() {
    vector<string> items;
    items.push_back("1 - x\'\' + omega^2 * (1 + a * cos(t)) * x = 0");
    items.push_back("2 - x\'\' + (k + delta * cos(t)) * x = 0");
    items.push_back("3 - Periodic perturbation of Euler-Poinsot case");
    items.push_back("4 - Periodic perturbation of Joukowskii-Volterra case");
    items.push_back("5 - Multi-threading options");
    items.push_back("0 - Exit");

    int selected = -1;
    while (selected != 0) {
        selected = -1;
        system("clear");
        for (int n = 0; n < (int)items.size(); n++) {
            printf("%s\n", items[n].c_str());
        }
        printf("Select: ");
        
        int res = scanf("%d", &selected);
        if (res != 1) {
            clearStream();
            continue;
        }
        clearStream();
        switch (selected) {
        case 1:

            break;
        case 2:

            break;
        case 3:

            break;
        case 4:

            break;
        case 5:
            multitreadingOptions();
            break;
        }

    }

}

void saveToFile(const char * fName, valarray<double> & x, valarray<double> & y, 
        valarray<valarray<double> > & T) {
    FILE *f = fopen(fName, "w");
    int N_x = (int)x.size();
    int N_y = (int)y.size();

    for (int i = 0; i < N_x; i++) {
        fprintf(f, "\t%.15lg", x[i]);
    }
    fprintf(f, "\n");

    for (int j = 0; j < N_y; j++) {
        fprintf(f, "%.15lg", y[j]);
        for (int i = 0; i < N_x; i++) {
            fprintf(f, "\t%.15lg", T[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void case1() {
    double omega_min = 0.0;
    double omega_max = 2.5;
    double e_min = 0.0;
    double e_max = 2.0;

    int N_omega = 1001;
    int N_e = 1001;
    
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
    double ic1[2] = { 1.0, 0.0 }, ic2[2] = { 0.0, 1.0 };
   
    int count = 0;
    
    double start = omp_get_wtime();
    #pragma omp parallel
    {
        RungeKutta method;
        Mathieu eqs;
        method.init(&eqs);

        double img1[2], img2[2];
        double A[2][2];
        #pragma omp for
        for (int i = 0; i < N_omega; i++) {
            for (int j = 0; j < N_e; j++) {
                eqs.setParameter(0, omega[i]);
                eqs.setParameter(1, e[j]);

                method.map(&eqs, 0.0, ic1, img1, N_t, T);
                method.map(&eqs, 0.0, ic2, img2, N_t, T);
            
                A[0][0] = img1[0]; A[0][1] = img2[0]; 
                A[1][0] = img1[1]; A[1][1] = img2[1];

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
}

void case2() {
    double omega_min = -0.5;
    double omega_max = 2.5;
    double e_min = 0.0;
    double e_max = 1.0;

    int N_omega = 1001;
    int N_e = 1001;
    
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
    double ic1[2] = { 1.0, 0.0 }, ic2[2] = { 0.0, 1.0 };
   
    int count = 0;
    
    double start = omp_get_wtime();
    #pragma omp parallel
    {
        RungeKutta method;
        Mathieu2 eqs;
        method.init(&eqs);

        double img1[2], img2[2];
        double A[2][2];
        #pragma omp for
        for (int i = 0; i < N_omega; i++) {
            for (int j = 0; j < N_e; j++) {
                eqs.setParameter(0, omega[i]);
                eqs.setParameter(1, e[j]);

                method.map(&eqs, 0.0, ic1, img1, N_t, T);
                method.map(&eqs, 0.0, ic2, img2, N_t, T);
            
                A[0][0] = img1[0]; A[0][1] = img2[0]; 
                A[1][0] = img1[1]; A[1][1] = img2[1];

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
}

void case3() {

    double omega_min = 0.0;
    double omega_max = 6.0;
    double dJ2_min = 0.0;
    double dJ2_max = 2.0;

    int N_omega = 501;
    int N_dJ2 = 501;
    
    double dOmega = (omega_max - omega_min) / (N_omega - 1);
    double ddJ2 = (dJ2_max - dJ2_min) / (N_dJ2 - 1);

    valarray<double> dJ2(N_dJ2);
    valarray<double> omega(N_omega);
    valarray<valarray<double> > tr(N_omega), det(N_omega);
    for (int i = 0; i < N_omega; i++) {
        tr[i].resize(N_dJ2, 0.0);
        det[i].resize(N_omega, 0.0);
    }

    for (int i = 0; i < N_omega; i++) {
        omega[i] = omega_min + i * dOmega;
    }

    for (int j = 0; j < N_dJ2; j++) {
        dJ2[j] = dJ2_min + j * ddJ2;
    }

    double T = 1.0;
    int N_t = 10000;
    printf("Integration Step: %lg\n", T / N_t);
    double ic1[2] = { 1.0, 0.0 }, ic2[2] = { 0.0, 1.0 };
   
    int count = 0;
    
    double start = omp_get_wtime();
    #pragma omp parallel
    {
        RungeKutta method;
        LiouvilleHill eqs;
        method.init(&eqs);
        eqs.setParameter(0, 2.0);
        eqs.setParameter(1, 3.0);
        eqs.setParameter(2, 4.0);

        double img1[2], img2[2];
        double A[2][2];
        #pragma omp for
        for (int i = 0; i < N_omega; i++) {
            for (int j = 0; j < N_dJ2; j++) {
                eqs.setParameter(3, dJ2[j]);
                eqs.setParameter(4, omega[i]);

                method.map(&eqs, 0.0, ic1, img1, N_t, T);
                method.map(&eqs, 0.0, ic2, img2, N_t, T);
            
                A[0][0] = img1[0]; A[0][1] = img2[0]; 
                A[1][0] = img1[1]; A[1][1] = img2[1];

                det[i][j] = A[0][0] * A[1][1] - A[1][0] * A[0][1];
                tr[i][j] = fabs(A[0][0] + A[1][1]);
            }
            
            #pragma omp critical
            {
                count += N_dJ2;
                printf("Ready %6d from %6d\n", count, N_omega * N_dJ2);
            }
        }
    }
    printf("Time: %lg\n", omp_get_wtime() - start);
    saveToFile("trace.csv", omega, dJ2, tr);
    saveToFile("det.csv", omega, dJ2, det);
}

void case4() {
    double OMEGA_N_CONST = 200.0 * M_PI;

    double omega_min = 0.02;
    double omega_max = 0.22;
    double dJ2_min = 0.0;
    double dJ2_max = 0.5;

    int N_omega = 16001;
    int N_dJ2 = 1001;
    
    double ddJ2 = (dJ2_max - dJ2_min) / (N_dJ2 - 1);

    valarray<double> dJ2(N_dJ2);
    valarray<double> omega(N_omega);
    valarray<valarray<double> > tr(N_omega), det(N_omega);
    for (int i = 0; i < N_omega; i++) {
        tr[i].resize(N_dJ2, 0.0);
        det[i].resize(N_dJ2, 0.0);
    }

    double p1 = log10(omega_min);
    double p2 = log10(omega_max);
    //double dOmega = (omega_max - omega_min) / (N_omega - 1);
    double dOmega = (p2 - p1) / (N_omega - 1);
    for (int i = 0; i < N_omega; i++) {
        omega[i] = pow(10, p1 + i * dOmega);
    }

    for (int j = 0; j < N_dJ2; j++) {
        dJ2[j] = dJ2_min + j * ddJ2;
    }

    double T = 1.0;
    //printf("Integration Step: %lg\n", T / N_t);
    double ic1[2] = { 1.0, 0.0 }, ic2[2] = { 0.0, 1.0 };
   
    int count = 0;
    
    double start = omp_get_wtime();
    #pragma omp parallel
    {
        RungeKutta method;
        LiouvilleHill2 eqs;
        method.init(&eqs);
        eqs.setParameter(0, 2.0);
        eqs.setParameter(1, 3.0);
        eqs.setParameter(2, 4.0);
        eqs.setParameter(5, 0.0);

        double img1[2], img2[2];
        double A[2][2];
        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < N_omega; i++) {
            for (int j = 0; j < N_dJ2; j++) {
                eqs.setParameter(3, dJ2[j]);
                eqs.setParameter(4, omega[i]);

                int N_t = (int)ceil(OMEGA_N_CONST / omega[i]);

                method.map(&eqs, 0.0, ic1, img1, N_t, T);
                method.map(&eqs, 0.0, ic2, img2, N_t, T);
            
                A[0][0] = img1[0]; A[0][1] = img2[0]; 
                A[1][0] = img1[1]; A[1][1] = img2[1];

                det[i][j] = A[0][0] * A[1][1] - A[1][0] * A[0][1];
                tr[i][j] = fabs(A[0][0] + A[1][1]);
            }
            
            #pragma omp critical
            {
                count += N_dJ2;
                printf("Ready %6d from %6d\n", count, N_omega * N_dJ2);
            }
        }
    }
    printf("Time: %lg\n", omp_get_wtime() - start);
    saveToFile("trace.csv", omega, dJ2, tr);
    saveToFile("det.csv", omega, dJ2, det);
}
