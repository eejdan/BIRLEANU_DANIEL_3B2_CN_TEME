#include <cmath> 
#include <math.h>
#include <stdio.h>

#include "mpu.h"
#include "tantester.h"

// VV ca sa taca vscode-ul
#ifndef M_PI
    #include <corecrt_math_defines.h>
#endif

class TanUtilBase {
protected:
    static double reduceToTanPeriod(double x) { // function body facut cu LLM
        double reduced = fmod(x, M_PI);
        
        if (reduced > M_PI_2) {
            reduced -= M_PI;
        } else if (reduced <= -M_PI_2) {
            reduced += M_PI;
        }
        
        return reduced;
    }
};
class MyTanLentz : public TanUtilBase {
public:
    static double myTan(double x) {
        double mpu = getMachinePrecisionUnit();
        if(x <= -M_PI_2 || x >= M_PI_2) {
            x = reduceToTanPeriod(x);
        }
        return myTanLentzImpl(x, mpu);
    }
    static double myTanGivenPrecision(double x, double epsilon) {
        if(x <= -M_PI_2 || x >= M_PI_2) {
            x = reduceToTanPeriod(x);
        }
        return myTanLentzImpl(x, epsilon);
    }
private:
    MyTanLentz() {}
    static double getScalarBi(int i) { return i == 0 ? 0.0 : (double)(2 * i - 1); }
    static double getScalarAi(int i, double x) { return i < 2 ? x : -pow(x, 2); }
    static double myTanLentzImpl(double x, double epsilon) {
        double f = getScalarBi(0);
        double prec = pow(10.0, -12);
        if(f == 0) f = prec;
        double C = f;
        double D = 0;
        int j = 1;
        double delta;
        do {
            D = getScalarBi(j) + getScalarAi(j, x) * D;
            if(D == 0) D = prec;
            C = getScalarBi(j) + getScalarAi(j, x) / C;
            if(C == 0) C = prec;
            D = 1 / D;
            delta = C * D;
            f *= delta;
            j++;
        } while(abs(delta - 1.0) >= epsilon);
        return f;
    }
};
class MyTanPolinomial : public TanUtilBase {
public:
    static double myTan(double x) {
        if(x <= -M_PI_2 || x >= M_PI_2) {
            x = reduceToTanPeriod(x);
        }
        if(x <= -M_PI_4) {
            return -myTanPolinomialImpl(-x);
        } else if(x >= M_PI_4) {
            return 1.0 / myTanPolinomialImpl(M_PI_2 - x);
        }
        return myTanPolinomialImpl(x);
    }
private:
    static volatile const double coefficients[5];
    static double myTanPolinomialImpl(double x) {
        double result = 0.0;
        for(int i = 0; i < 5; i++) {
            result += coefficients[i] * pow(x, (i+1)*2-1);
        }
        return result;
    }
};
volatile const double MyTanPolinomial::coefficients[5] = {1.0, 1.0/3.0, 2.0/15.0, 17.0/315.0, 62.0/2835.0};

int main() {
    double testValue = 77.0 * M_PI / 6.0;
    printf("tan(pi/6) = \t\t\t\t%.21f\n", tan(testValue));
    printf("myTanLentz(pi/6) = \t\t\t%.21f\n", MyTanLentz::myTan(testValue));
    printf("myTanPolinomial(pi/6) = \t\t%.21f\n", MyTanPolinomial::myTan(testValue));
    printf("| tan(pi/6) - myTanLentz(pi/6) | =\t %.21f\n", abs(tan(testValue) - MyTanLentz::myTan(testValue)));
    printf("| tan(pi/6) - myTanPolinomial(pi/6) | =  %.21f\n", abs(tan(testValue) - MyTanPolinomial::myTan(testValue)));

    printf("\n");

    TanTestResult lentzResult = testTanApproximation(MyTanLentz::myTan);
    printTanTestResult(lentzResult, "Lentz Continued Fraction");

    TanTestResult polyResult = testTanApproximation(MyTanPolinomial::myTan);
    printTanTestResult(polyResult, "Polynomial (Taylor)");

    return 0;
}