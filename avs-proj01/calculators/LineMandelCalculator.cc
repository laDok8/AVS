/**
 * @file LineMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>


#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
    data = (int *)(malloc(height * width * sizeof(int)));
    realLine = (float *)(malloc(width * sizeof(float)));
    imagLine = (float *)(malloc(width * sizeof(float)));
    realLineStart = (float *)(malloc(width * sizeof(float)));

    for (int j = 0; j < width; j++) {
        const float x = x_start + j * dx; // current real value
        realLineStart[j] = x;
    }

}

LineMandelCalculator::~LineMandelCalculator() {
    free(data);
    data = NULL;
    free(realLine);
    realLine = NULL;
    free(imagLine);
    imagLine = NULL;
    free(realLineStart);
    realLineStart = NULL;
}

template <typename T>
static inline void mandelbrotLine(int width, int *pdata, T *realLine, T *imagLine, T *realLineStart, T imagLineStart)
{
    #pragma omp simd linear(p:1) uniform(imagLineStart) private(r2,i2)
    for (int p= 0; p < width; p++)
    {
        T r2 = realLine[p] * realLine[p];
        T i2 = imagLine[p] * imagLine[p];

        if (r2 + i2 > 4.0f) {
            continue;
        }

        imagLine[p] = 2.0f * realLine[p] * imagLine[p] + imagLineStart;
        realLine[p] = r2 - i2 + realLineStart[p];
        pdata[p]++;
    }

}



int * LineMandelCalculator::calculateMandelbrot () {
    int *pdata = data;
    for (int i = 0; i < height; i++)
    {
        float y = y_start + i * dy; // current imaginary value

        //init assign
        #pragma omp simd linear(j:1) uniform(y)
        for (int j = 0; j < width; j++)
        {
            imagLine[j] = y;
            float x = x_start + j * dx; // current real value
            realLine[j] = realLineStart[j] = x;
        }

        for (int i = 0; i < limit; ++i)
        {
            mandelbrotLine(width,pdata,realLine,imagLine,realLineStart,y);
        }
        pdata+=width;
    }
    return data;
}
