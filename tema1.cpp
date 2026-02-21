#include <cmath>
#include <stdio.h>
#include <stdexcept>
// compilat fara argumente (">g++ tema1.cpp"); folosesc gcc version 6.3.0 (MinGW.org GCC-6.3.0-1)
// Ex 1.
double getMachinePrecisionUnit() {
    int m=1;
     
    while(m < 1000) {
        double x = 1.0 + pow(10.0, -m);
        if (x == 1.0) return pow(10.0, -(m-1));
        m++;
    }

    throw std::runtime_error("Nu s-a gasit precizia masina pana in limita");
}

// Ex 2.
void canAdditionAssociativityBreak() {
    double mpu = getMachinePrecisionUnit();

    double x = 1.0;
    volatile double y = mpu/10.0;
    volatile double z = mpu/10.0;
    printf("(x, y, z) = (%e, %e, %e)\n", x, y, z);

    printf("testez (x+y)+z == x+(y+z) => %s\n", ((x + y) + z) == (x  + (y + z)) ? "este asociativ" : "nu e asociativ");
}
void canMultiplicationAssociativityBreak() {
    double x = ((double) rand() / RAND_MAX);
    double y = ((double) rand() / RAND_MAX);
    double z = ((double) rand() / RAND_MAX);
    while(x * (y * z) == (x * y) * z) {
        y = ((double) rand() / RAND_MAX);
        z = ((double) rand() / RAND_MAX);
    }
    printf("(x, y, z) = (%e, %e, %e)\n", x, y, z);
    printf("testez (x*y)*z == x*(y*z) => %s\n", ((x * y) * z) == (x  * (y * z)) ? "este asociativ" : "nu e asociativ");
    // ar trebui sa gaseasca in sub o secunda, daca nu, rerun
}
// Ex 3.


int main() {
    double mpu = getMachinePrecisionUnit();
    printf("Machine precision unit: %e\n", mpu);
    canAdditionAssociativityBreak();
    canMultiplicationAssociativityBreak();
    return 0;
}