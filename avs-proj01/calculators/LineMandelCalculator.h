/**
 * @file LineMandelCalculator.h
 * @author Ladislav Dokoupil (xdokou14@fit.vutbr.cz)
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */

#include <BaseMandelCalculator.h>

class LineMandelCalculator : public BaseMandelCalculator
{
public:
    LineMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~LineMandelCalculator();
    int *calculateMandelbrot();

private:
    int *data;
    float *realLine,*imagLine;
    float *realLineStart;

};