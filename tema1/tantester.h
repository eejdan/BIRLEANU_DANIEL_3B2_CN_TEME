#ifndef TANTESTER_H
#define TANTESTER_H

#include <functional>
// fisier generat cu LLM
struct TanTestResult {
    double maxAbsError;
    double avgAbsError;
    double maxRelError;
    double avgRelError;
    double approxTimeMs;
    double stdTanTimeMs;
    int numTests;
};

TanTestResult testTanApproximation(
    std::function<double(double)> approxTan,
    int numTests = 10000,
    double rangeMin = -1.5,
    double rangeMax = 1.5
);

void printTanTestResult(const TanTestResult& result, const char* name);

#endif
