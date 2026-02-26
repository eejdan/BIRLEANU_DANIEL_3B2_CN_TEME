#include "tantester.h"

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>

//fisier generat cu LLM

TanTestResult testTanApproximation(
    std::function<double(double)> approxTan,
    int numTests,
    double rangeMin,
    double rangeMax
) {
    TanTestResult result = {};
    result.numTests = numTests;
    result.maxAbsError = 0.0;
    result.avgAbsError = 0.0;
    result.maxRelError = 0.0;
    result.avgRelError = 0.0;

    // Generate test values evenly spaced in [rangeMin, rangeMax]
    double step = (rangeMax - rangeMin) / (numTests - 1);

    // --- Time the standard tan ---
    auto startStd = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < numTests; i++) {
        double x = rangeMin + i * step;
        volatile double y = std::tan(x); // volatile to prevent optimizing away
        (void)y;
    }
    auto endStd = std::chrono::high_resolution_clock::now();
    result.stdTanTimeMs = std::chrono::duration<double, std::milli>(endStd - startStd).count();

    // --- Time the approximation and compute errors ---
    auto startApprox = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < numTests; i++) {
        double x = rangeMin + i * step;
        volatile double y = approxTan(x);
        (void)y;
    }
    auto endApprox = std::chrono::high_resolution_clock::now();
    result.approxTimeMs = std::chrono::duration<double, std::milli>(endApprox - startApprox).count();

    // --- Compute errors in a separate pass (so timing is clean) ---
    for (int i = 0; i < numTests; i++) {
        double x = rangeMin + i * step;
        double expected = std::tan(x);
        double approx = approxTan(x);

        double absErr = std::fabs(expected - approx);
        double relErr = (expected != 0.0) ? std::fabs(absErr / expected) : absErr;

        result.avgAbsError += absErr;
        result.avgRelError += relErr;

        if (absErr > result.maxAbsError) result.maxAbsError = absErr;
        if (relErr > result.maxRelError) result.maxRelError = relErr;
    }

    result.avgAbsError /= numTests;
    result.avgRelError /= numTests;

    return result;
}

void printTanTestResult(const TanTestResult& result, const char* name) {
    printf("=== Tan Approximation Test: %s ===\n", name);
    printf("  Test points:       %d\n", result.numTests);
    printf("  Max abs error:     %.15e\n", result.maxAbsError);
    printf("  Avg abs error:     %.15e\n", result.avgAbsError);
    printf("  Max rel error:     %.15e\n", result.maxRelError);
    printf("  Avg rel error:     %.15e\n", result.avgRelError);
    printf("  Approx time:       %.4f ms\n", result.approxTimeMs);
    printf("  std::tan time:     %.4f ms\n", result.stdTanTimeMs);
    printf("  Speedup (std/approx): %.2fx\n", result.stdTanTimeMs / result.approxTimeMs);
    printf("==========================================\n\n");
}
