/**
 * @file BatchMandelCalculator.h
 * @author LADISLAV DOKOUPIL <xdokou14@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 6.11.2022
 */
#ifndef BATCHMANDELCALCULATOR_H
#define BATCHMANDELCALCULATOR_H

#include <BaseMandelCalculator.h>

class BatchMandelCalculator : public BaseMandelCalculator
{
public:
    BatchMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~BatchMandelCalculator();
    int * calculateMandelbrot();

private:
    int *data;
    float *realBlock, *imagBlock, *realBlockStart;
    const int batch_size = 256;
};

#endif