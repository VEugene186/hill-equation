#include "Equation.h"

Equation::Equation(int dim, int numberOfParameters) : 
        dim_(dim), numberOfParameters_(numberOfParameters), parameters_(nullptr) {
    if (numberOfParameters_ > 0) {
        parameters_ = new double[numberOfParameters_];
    }
}

Equation::~Equation() {
    if (parameters_ != nullptr) {
        delete [] parameters_;
    }
}

int Equation::getDim() const {
    return dim_;
}

int Equation::getNumberOfParameters() const {
    return numberOfParameters_;
}

double Equation::getParameter(int n) const {
    return parameters_[n];
}

void Equation::setParameter(int n, double value) {
    parameters_[n] = value;
}
