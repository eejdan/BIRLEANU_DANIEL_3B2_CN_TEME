#include <stdio.h>

#include "mpu.h"

int main() {
    double mpu = getMachinePrecisionUnit();
    printf("Precizia masina: %e\n", mpu);
    return 0;
}