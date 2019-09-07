#include "RungeKutta.h"
#include <cstdlib>
#include <cmath>

RungeKutta::RungeKutta() : dim_(0), 
        k1_(nullptr), k2_(nullptr), k3_(nullptr), k4_(nullptr), q_(nullptr), qTmp_(nullptr) {

}

RungeKutta::~RungeKutta() {
    if (k1_ != nullptr) {
        delete [] k1_;
    }
    if (k2_ != nullptr) {
        delete [] k2_;
    }
    if (k3_ != nullptr) {
        delete [] k3_;
    }
    if (k4_ != nullptr) {
        delete [] k4_;
    }
    if (q_ != nullptr) {
        delete [] q_;
    }
    if (qTmp_ != nullptr) {
        delete [] qTmp_;
    }
}

void RungeKutta::init(const Equation *eq) {
    if (dim_ == eq->getDim()) return;
    dim_ = eq->getDim();
    if (dim_ > 0) {
        k1_ = new double[dim_];
        k2_ = new double[dim_];
        k3_ = new double[dim_];
        k4_ = new double[dim_];
        q_ = new double[dim_];
        qTmp_ = new double[dim_];
    }
}

void RungeKutta::makeStep(const Equation *eq, double t0, const double *q0, double *q1, double dt) {
    double dt_2 = 0.5 * dt;
    double dt_6 = dt / 6.0;
   
    eq->RHS(t0, q0, k1_);
    for (int i = 0; i < dim_; i++) {
        q_[i] = q0[i] + dt_2 * k1_[i]; 
    }
    
    eq->RHS(t0 + dt_2, q_, k2_);
    for (int i = 0; i < dim_; i++) {
        q_[i] = q0[i] + dt_2 * k2_[i]; 
    }

    eq->RHS(t0 + dt_2, q_, k3_);
    for (int i = 0; i < dim_; i++) {
        q_[i] = q0[i] + dt * k3_[i]; 
    }

    eq->RHS(t0 + dt, q_, k4_);
    for (int i = 0; i < dim_; i++) {
        q1[i] = q0[i] + dt_6 * (k1_[i] + 2.0 * (k2_[i] + k3_[i]) + k4_[i]); 
    }
}

void RungeKutta::map(const Equation *eq, double t0, const double *q0, double *q1, int N_steps, double T) {
    for (int i = 0; i < dim_; i++) {
        qTmp_[i] = q0[i];
    }
    double dt = T / N_steps;
    for (int n = 0; n < N_steps; n++) {
        makeStep(eq, t0, qTmp_, q1, dt);
        t0 += dt;
        for (int i = 0; i < dim_; i++) {
            qTmp_[i] = q1[i];
        }
    }
    //if (fabs(t0 - T) > 1e-10) exit(1);
}

