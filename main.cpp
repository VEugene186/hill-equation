#include <cstdio>
#include <cmath>
#include <valarray>
#include <vector>
#include <string>
#include "Mathieu.h"
#include "Mathieu2.h"
#include "LiouvilleHill.h"
#include "EulerPoisson.h"
#include "RungeKutta.h"
#include "EulerPoissonFixedPoints.h"
#include "EulerPoissonTracker.h"
#include <omp.h>

using namespace std;

void saveToFile(const char * fName, valarray<double> & x, valarray<double> & y, 
        valarray<valarray<double> > & T);

void startMenuLoop();

void multitreadingOptions();
void clearStream();

void MathieuEquation_Arnold();
void MathieuEquation_Kovacic();
void EulerPoinsotStability();
void JoukowskiiVolterraStability_Flow(); 
void JoukowskiiVolterraStability_Map(); 
void JoukowskiiVolterraStability_findFP(); 
void JoukowskiiVolterraStability_trackFP(); 
void JoukowskiiVolterraStability_trackMap(); 

int main(int argc, char * argv[]) {
    startMenuLoop();
    
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
    items.push_back("2 - x\'\' + (delta + epsilon * cos(t)) * x = 0");
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
            MathieuEquation_Arnold();
            break;
        case 2:
            MathieuEquation_Kovacic();
            break;
        case 3:
            EulerPoinsotStability();
            break;
        case 4:
            //JoukowskiiVolterraStability_Flow();
            //JoukowskiiVolterraStability_Map();
            //JoukowskiiVolterraStability_findFP();
            //JoukowskiiVolterraStability_trackFP();
            JoukowskiiVolterraStability_trackMap(); 
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

void inputMesh(const char * par1Name, const char * par2Name, int * par1N, double * par1min, double * par1max, 
                                                             int * par2N, double * par2min, double * par2max) {
    printf("%s, nodes: ", par1Name);
    scanf("%d", par1N);
    printf("%s, min  : ", par1Name);
    scanf("%lf", par1min);
    printf("%s, max  : ", par1Name);
    scanf("%lf", par1max);

    printf("%s, nodes: ", par2Name);
    scanf("%d", par2N);
    printf("%s, min  : ", par2Name);
    scanf("%lf", par2min);
    printf("%s, max  : ", par2Name);
    scanf("%lf", par2max);
}

void createMesh(int N, double left, double right, valarray<double> & mesh, bool logScale = false) {
    double step;
    mesh.resize(N);
    if (logScale) {
        double p1 = log10(left);
        double p2 = log10(right);
        step = (p2 - p1) / (N - 1);
        for (int i = 0; i < N; i++)
            mesh[i] = pow(10, p1 + i * step);
    }
    else {
        step = (right - left) / (N - 1);
        for (int i = 0; i < N; i++)
            mesh[i] = left + i * step;
    }
}

void allocate2d(int N1, int N2, valarray<valarray<double> > & A) {
    A.resize(N1);
    for (int i = 0; i < N1; i++) {
        A[i].resize(N2, 0.0);
    }
}

void MathieuEquation_Arnold() {
    double omega_min, omega_max;
    double a_min, a_max;
    int N_omega, N_a;

    inputMesh("omega", "a", &N_omega, &omega_min, &omega_max, &N_a, &a_min, &a_max);

    valarray<double> omega, a;
    createMesh(N_omega, omega_min, omega_max, omega);
    createMesh(N_a, a_min, a_max, a);

    valarray<valarray<double> > tr, det;
    allocate2d(N_omega, N_a, tr);
    allocate2d(N_omega, N_a, det);


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
            for (int j = 0; j < N_a; j++) {
                eqs.setParameter(0, omega[i]);
                eqs.setParameter(1, a[j]);

                method.map(&eqs, 0.0, ic1, img1, N_t, T);
                method.map(&eqs, 0.0, ic2, img2, N_t, T);
            
                A[0][0] = img1[0]; A[0][1] = img2[0]; 
                A[1][0] = img1[1]; A[1][1] = img2[1];

                det[i][j] = A[0][0] * A[1][1] - A[1][0] * A[0][1];
                tr[i][j] = fabs(A[0][0] + A[1][1]);
            }
            
            #pragma omp critical
            {
                count += N_a;
                printf("Ready %6d from %6d\n", count, N_omega * N_a);
            }
        }
    }
    printf("Time: %lg\n", omp_get_wtime() - start);
    saveToFile("trace.csv", omega, a, tr);
    saveToFile("det.csv", omega, a, det);
}

void MathieuEquation_Kovacic() {
    double delta_min, delta_max;
    double epsilon_min, epsilon_max;
    int N_delta, N_epsilon;
    
    inputMesh("delta", "epsilon", &N_delta, &delta_min, &delta_max, &N_epsilon, &epsilon_min, &epsilon_max);

    valarray<double> delta, epsilon;
    createMesh(N_delta, delta_min, delta_max, delta);
    createMesh(N_epsilon, epsilon_min, epsilon_max, epsilon);

    valarray<valarray<double> > tr, det;
    allocate2d(N_delta, N_epsilon, tr);
    allocate2d(N_delta, N_epsilon, det);

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
        for (int i = 0; i < N_delta; i++) {
            for (int j = 0; j < N_epsilon; j++) {
                eqs.setParameter(0, delta[i]);
                eqs.setParameter(1, epsilon[j]);

                method.map(&eqs, 0.0, ic1, img1, N_t, T);
                method.map(&eqs, 0.0, ic2, img2, N_t, T);
            
                A[0][0] = img1[0]; A[0][1] = img2[0]; 
                A[1][0] = img1[1]; A[1][1] = img2[1];

                det[i][j] = A[0][0] * A[1][1] - A[1][0] * A[0][1];
                tr[i][j] = fabs(A[0][0] + A[1][1]);
            }
            
            #pragma omp critical
            {
                count += N_epsilon;
                printf("Ready %6d from %6d\n", count, N_delta * N_epsilon);
            }
        }
    }
    printf("Time: %lg\n", omp_get_wtime() - start);
    saveToFile("trace.csv", delta, epsilon, tr);
    saveToFile("det.csv", delta, epsilon, det);
}



void EulerPoinsotStability() {
    double OMEGA_N_CONST = 200.0 * M_PI;

    double Omega_min, Omega_max;
    double dJ2_min, dJ2_max;
    int N_Omega, N_dJ2;

    inputMesh("Omega", "dJ2", &N_Omega, &Omega_min, &Omega_max, &N_dJ2, &dJ2_min, &dJ2_max);

    valarray<double> Omega, dJ2;
    createMesh(N_Omega, Omega_min, Omega_max, Omega, true);
    createMesh(N_dJ2, dJ2_min, dJ2_max, dJ2);

    valarray<valarray<double> > tr, det;
    allocate2d(N_Omega, N_dJ2, tr);
    allocate2d(N_Omega, N_dJ2, det);

    double T = 1.0;
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
        eqs.setParameter(5, 0.0);

        double img1[2], img2[2];
        double A[2][2];
        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < N_Omega; i++) {
            for (int j = 0; j < N_dJ2; j++) {
                eqs.setParameter(3, dJ2[j]);
                eqs.setParameter(4, Omega[i]);

                int N_t = (int)ceil(OMEGA_N_CONST / Omega[i]);

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
                printf("Ready %6d from %6d\n", count, N_Omega * N_dJ2);
            }
        }
    }
    printf("Time: %lg\n", omp_get_wtime() - start);
    saveToFile("trace.csv", Omega, dJ2, tr);
    saveToFile("det.csv", Omega, dJ2, det);
}

void JoukowskiiVolterraStability_Flow() {
    double M0[3] = { 0.536357900134574, 0.589984383362725, 3.37108864269863 };//1 stable point
    //double M0[3] = { 0.5363579, 0.589984383362725, 3.37108864269863 };//1 stable point perturbed
    //double M0[3] = { 0.792282734739558, 2.94208261906373, -1.64816198565656 };//2 unstable point
    //double M0[3] = { 0.79228273473, 2.94208261906373, -1.64816198565656 };//2 unstable point perturbed
    EulerPoisson eqs;
    eqs.setParameter(0, 2.0);
    eqs.setParameter(1, 3.0);
    eqs.setParameter(2, 4.0);
    eqs.setParameter(3, 0.3);
    eqs.setParameter(4, 0.2);
    eqs.setParameter(5, 0.4);
    eqs.setParameter(6, 2.0 * M_PI);
    eqs.setParameter(7, 0.01);
    RungeKutta method;
    method.init(&eqs);

    double T = 200.0;
    int N = 20001;
    double dt = T / (N - 1);
    valarray<double> t(N);
    valarray<double*> sol(N);
    for (int i = 0; i < N; i++) {
        sol[i] = new double[3];
    }
    t[0] = 0.0;
    sol[0][0] = M0[0];
    sol[0][1] = M0[1];
    sol[0][2] = M0[2];

    for (int i = 0; i < N - 1; i++) {
        t[i + 1] = t[i] + dt;
        method.makeStep(&eqs, t[i], sol[i], sol[i + 1], dt);
    }

    FILE * f = fopen("solution.csv", "w");
    for (int i = 0; i < N; i++) {
        fprintf(f, "%.15lg\t%.15lg\t%.15lg\t%.15lg\n", t[i], sol[i][0], sol[i][1], sol[i][2]);
    }
    fclose(f);        
    for (int i = 0; i < N; i++) {
        delete [] sol[i];
    }
}

void JoukowskiiVolterraStability_Map() {
    //double M0[3] = { 0.536357900134574, 0.589984383362725, 3.37108864269863 };//1 stable point
    //double M0[3] = { 0.5363579, 0.589984383362725, 3.37108864269863 };//1 stable point perturbed
    //double M0[3] = { 0.792282734739558, 2.94208261906373, -1.64816198565656 };//2 unstable point
    double M0[3] = { 0.79228273473, 2.94208261906373, -1.64816198565656 };//2 unstable point perturbed
    EulerPoisson eqs;
    eqs.setParameter(0, 2.0);
    eqs.setParameter(1, 3.0);
    eqs.setParameter(2, 4.0);
    eqs.setParameter(3, 0.3);
    eqs.setParameter(4, 0.2);
    eqs.setParameter(5, 0.4);
    eqs.setParameter(6, 2.0 * M_PI);
    eqs.setParameter(7, 0.0);
    RungeKutta method;
    method.init(&eqs);

    double T = 1.0;
    int N = 20001;
    valarray<double> t(N);
    valarray<double*> sol(N);
    for (int i = 0; i < N; i++) {
        sol[i] = new double[3];
    }
    t[0] = 0.0;
    sol[0][0] = M0[0];
    sol[0][1] = M0[1];
    sol[0][2] = M0[2];

    for (int i = 0; i < N - 1; i++) {
        t[i + 1] = t[i] + T;
        method.map(&eqs, t[i], sol[i], sol[i + 1], 100, T);
    }

    FILE * f = fopen("solution.csv", "w");
    for (int i = 0; i < N; i++) {
        fprintf(f, "%.15lg\t%.15lg\t%.15lg\t%.15lg\n", t[i], sol[i][0], sol[i][1], sol[i][2]);
    }
    fclose(f);        
    for (int i = 0; i < N; i++) {
        delete [] sol[i];
    }
}

void JoukowskiiVolterraStability_findFP() {
    //double M0[3] = { 0.536357900134574, 0.589984383362725, 3.37108864269863 };//1 stable point
    //double M0[3] = { 0.5363579, 0.589984383362725, 3.37108864269863 };//1 stable point perturbed
    //double M0[3] = { 0.792282734739558, 2.94208261906373, -1.64816198565656 };//2 unstable point
    //double M0[3] = { 0.79228273473, 2.94208261906373, -1.64816198565656 };//2 unstable point perturbed
    double M0[3] = { -3.43273007941469, -0.316911509971415, -0.340486966577532 };//4 stable point
    
    double M1[3], phi1, z1;

    double g = sqrt(M0[0] * M0[0] + M0[1] * M0[1] + M0[2] * M0[2]);

    EulerPoisson eqs;
    eqs.setParameter(0, 2.0);
    eqs.setParameter(1, 3.0);
    eqs.setParameter(2, 4.0);
    eqs.setParameter(3, 0.3);
    eqs.setParameter(4, 0.2);
    eqs.setParameter(5, 0.4);
    eqs.setParameter(6, 0.5);
    eqs.setParameter(7, 0.2);

    EulerPoissonFixedPoints epfp;
    epfp.find(&eqs, M0, M1, &phi1, &z1);

    printf("%.15lg\t%.15lg\t%.15lg\n\n", M0[0], M0[1], M0[2]);
    printf("%.15lg\t%.15lg\t%.15lg\n", M1[0], M1[1], M1[2]);
    printf("%.15lg\t%.15lg\n", atan2(M1[1], M1[0]), M1[2] / g);
    getchar();
}

void  JoukowskiiVolterraStability_trackFP() {
    //double M0[3] = { 0.536357900134574, 0.589984383362725, 3.37108864269863 };//1 stable point
    //double M0[3] = { 0.5363579, 0.589984383362725, 3.37108864269863 };//1 stable point perturbed
    //double M0[3] = { 0.792282734739558, 2.94208261906373, -1.64816198565656 };//2 unstable point
    //double M0[3] = { 0.79228273473, 2.94208261906373, -1.64816198565656 };//2 unstable point perturbed
    double M0[3] = { -3.43273007941469, -0.316911509971415, -0.340486966577532 };//4 stable point

    double dJ2_min = 0.0;
    double dJ2_max = 1.0;
    int N = 501;
    valarray<double*> points;
    valarray<double> dJ2, det(N), tr(N);
    
    EulerPoisson eqs;
    eqs.setParameter(0, 2.0);
    eqs.setParameter(1, 3.0);
    eqs.setParameter(2, 4.0);
    eqs.setParameter(3, 0.3);
    eqs.setParameter(4, 0.2);
    eqs.setParameter(5, 0.4);
    eqs.setParameter(6, 0.5);
    eqs.setParameter(7, 0.0);

    EulerPoissonTracker ept;
    ept.track(&eqs, M0, dJ2_min, dJ2_max, N, points, dJ2, tr, det);

    for (int i = 0; i < N; i++) {
        printf("%8lg | ", dJ2[i]);
        printf("%.15lg\t%.15lg\t%.15lg\t%.15lg\t%.15lg\n", 
                points[i][0], points[i][1], points[i][2], 
                points[i][3], points[i][4]);
    }
    getchar();

    
    for (int i = 0; i < N; i++) {
        delete [] points[i];
    }


}

void JoukowskiiVolterraStability_trackMap() {
    //double M0[3] = { 0.536357900134574, 0.589984383362725, 3.37108864269863 };//1 stable point
    //double M0[3] = { 0.5363579, 0.589984383362725, 3.37108864269863 };//1 stable point perturbed
    //double M0[3] = { 0.792282734739558, 2.94208261906373, -1.64816198565656 };//2 unstable point
    //double M0[3] = { 0.79228273473, 2.94208261906373, -1.64816198565656 };//2 unstable point perturbed
    double M0[3] = { -3.43273007941469, -0.316911509971415, -0.340486966577532 };//4 stable point

    double Omega_min, Omega_max;
    double dJ2_min, dJ2_max;
    int N_Omega, N_dJ2;

    inputMesh("Omega", "dJ2", &N_Omega, &Omega_min, &Omega_max, &N_dJ2, &dJ2_min, &dJ2_max);

    valarray<double> Omega, dJ2;
    createMesh(N_Omega, Omega_min, Omega_max, Omega, true);
    createMesh(N_dJ2, dJ2_min, dJ2_max, dJ2);

    valarray<valarray<double*> > points(N_Omega);
    valarray<valarray<double> > tr;
    valarray<valarray<double> > det;
    allocate2d(N_Omega, N_dJ2, tr);
    allocate2d(N_Omega, N_dJ2, det);

    int count = 0;
    double start = omp_get_wtime();
    #pragma omp parallel
    {
        EulerPoissonTracker ept;
        EulerPoisson eqs;
        eqs.setParameter(0, 2.0);
        eqs.setParameter(1, 3.0);
        eqs.setParameter(2, 4.0);
        eqs.setParameter(3, 0.3);
        eqs.setParameter(4, 0.2);
        eqs.setParameter(5, 0.4);
        //eqs.setParameter(6, 0.5);
        //eqs.setParameter(7, 0.0);
        valarray<double> tmp_dJ2;

        #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < N_Omega; i++) {
            eqs.setParameter(6, Omega[i]);
            ept.track(&eqs, M0, dJ2_min, dJ2_max, N_dJ2, points[i], tmp_dJ2, tr[i], det[i]);

            
            #pragma omp critical
            {
                count += N_dJ2;
                printf("Ready %6d from %6d\n", count, N_Omega * N_dJ2);
            }
        }
    }
    printf("Time: %lg\n", omp_get_wtime() - start);
    

    saveToFile("trace.csv", Omega, dJ2, tr);
    saveToFile("det.csv", Omega, dJ2, det);

}
