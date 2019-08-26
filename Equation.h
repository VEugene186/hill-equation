#ifndef EQUATION_H
#define EQUATION_H

class Equation {
public:
    Equation(int dim, int numberOfParameters);
    virtual ~Equation();

    virtual void RHS(double t, const double * q, double * dq) const = 0;

    int getDim() const;
    int getNumberOfParameters() const;
    double getParameter(int n) const;
    void setParameter(int n, double value);
protected:
    int dim_;
    int numberOfParameters_;
    double * parameters_;
};

#endif
