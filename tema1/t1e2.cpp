#include <cmath>
#include <stdio.h>

#include "mpu.h"

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


int main() {
    canAdditionAssociativityBreak();
    canMultiplicationAssociativityBreak();
    return 0;
}