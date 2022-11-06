/**
 * @file LineMandelCalculator.cc
 * @author LADISLAV DOKOUPIL <xdokou14@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 6.11.2022
 */

#include <algorithm>
#include <immintrin.h>

#include "LineMandelCalculator.h"

// malloc templates taken from school lab exercise
template<class T>
T* allocateMemory(size_t size)
{
    return ((T *) _mm_malloc(size * sizeof(T), 64));
}
template<class T>
void freeMemory(T* array)
{
    _mm_free(array);
    array = NULL;
}

LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{
    data = allocateMemory<int>( height * width * sizeof(int));
    realLine = allocateMemory<float>(width * sizeof(float));
    imagLine = allocateMemory<float>(width * sizeof(float));
    realLineStart = allocateMemory<float>(width * sizeof(float));
}

LineMandelCalculator::~LineMandelCalculator() {
    freeMemory(data);
    freeMemory(realLine);
    freeMemory(imagLine);
    freeMemory(realLineStart);
    }

int * LineMandelCalculator::calculateMandelbrot () {
    int *pdata = data;
    float *prealLine = realLine;
    float *pimagLine = imagLine;
    float *prealLineStart = realLineStart;

    int halfHeight = height / 2;
    for (int i = 0; i < halfHeight; i++)
    {
        float y = y_start + i * dy; // current imaginary value

        //init assign
        #pragma omp simd aligned(prealLine,pimagLine,prealLineStart)
        for (int j = 0; j < width; j++)
        {
            imagLine[j] = y;
            float x = x_start + j * dx; // current real value
            realLine[j] = realLineStart[j] = x;
        }

        for (int li = 0; li < limit; ++li)
        {
            //line
            int doneCount = 0;
            #pragma omp simd aligned(pdata,prealLine,pimagLine,prealLineStart) reduction(+:doneCount)
            for (int p= 0; p < width; p++)
            {
                float r2 = prealLine[p] * prealLine[p];
                float i2 = pimagLine[p] * pimagLine[p];

                if (r2 + i2 > 4.0f) {
                    doneCount++;
                    continue;
                }


                pimagLine[p] = 2.0f * prealLine[p] * pimagLine[p] + y;
                prealLine[p] = r2 - i2 + prealLineStart[p];
                pdata[p]++;
            }
            if(doneCount == width)
                break;
        }
        pdata+=width;
    }

    //data is symmetric
    for (int i = 0; i < halfHeight; i++)
    {
        for(int j = 0; j < width; j++)
        {
            data[(height -i - 1) * width + j] = data[i*width + j];
        }
    }

    return data;
}
