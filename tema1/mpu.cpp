#include <cmath>
#include <stdexcept>

double getMachinePrecisionUnit() {
    int m=1;
     
    while(m < 1000) {
        double x = 1.0 + pow(10.0, -m);
        if (x == 1.0) return pow(10.0, -(m-1));
        m++;
    }

    throw std::runtime_error("Nu s-a gasit precizia masina pana in limita");
}