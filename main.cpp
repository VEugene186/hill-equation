#include <cstdio>
#include <cmath>
#include <valarray>

using namespace std;

//Mathieu equation
double A(double t, double omega, double e) {
    return pow(omega, 2) * (1.0 + e * cos(t));
}

int main(int argc, char * argv[]) {
    double omega_min = 0.1;
    double omega_max = 2.0;
    double e_min = 0.0;
    double e_max = 1.0;

    int N_omega = 11;
    int N_e = 11;
    
    double dOmega = (omega_max - omega_min) / (N_omega - 1);
    double de = (e_max - e_min) / (N_e - 1);

    valarray<double> e(N_e);
    valarray<double> omega(N_omega);
    valarray<valarray<double> > tr(N_omega);
    for (int i = 0; i < N_omega; i++) {
        tr[i].resize(N_e, 0.0);
    }

    double T = 2.0 * M_PI;
    int N_t = 6000;
    double dt = T / N_t;
    printf("Integration Step: %lg\n", dt);
    for (int i = 0; i < N_omega; i++) {
        for (int j = 0; j < N_e; j++) {
            //loop on period
        }
    }
    
    return 0;
}
