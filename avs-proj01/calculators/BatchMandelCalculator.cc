/**
 * @file BatchMandelCalculator.cc
 * @author LADISLAV DOKOUPIL <xdokou14@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 6.11.2022
 */

#include <algorithm>
#include <stdlib.h>
#include <immintrin.h>

#include "BatchMandelCalculator.h"

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

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
    data = allocateMemory<int>( height * width * sizeof(int));
    realBlock = allocateMemory<float>(batch_size * sizeof(float));
    imagBlock = allocateMemory<float>(batch_size * sizeof(float));
    realBlockStart = allocateMemory<float>(batch_size * sizeof(float));
}

BatchMandelCalculator::~BatchMandelCalculator() {
    freeMemory(data);
    freeMemory(realBlock);
    freeMemory(imagBlock);
    freeMemory(realBlockStart);
}


int * BatchMandelCalculator::calculateMandelbrot () {
    int *pdata = data;
    float *prealBlock = realBlock;
    float *pimagBlock = imagBlock;
    float *prealBlockStart = realBlockStart;

    const int BLOCK = batch_size;
    const int halfHeight = height / 2;
    const int yBlocks = halfHeight / BLOCK;
    const int xBlocks = width / BLOCK;

    //y block
    for (int ty = 0; ty < yBlocks; ty++) {
        //x block
        for (int tx = 0; tx < xBlocks; tx++) {

            for (int y = 0; y < BLOCK; y++) {
                int effectiveY = ty * BLOCK + y;
                int dataBase = effectiveY * width;
                int *ppdata = pdata + dataBase + tx * BLOCK;

                //init values
                float startY = y_start + effectiveY * dy; // current imaginary value
                #pragma omp simd aligned(prealBlock,pimagBlock,prealBlockStart)
                for (int x = 0; x < BLOCK; x++) {
                    pimagBlock[x] = startY;
                    int effectiveX = tx * BLOCK + x;
                    float startX = x_start + effectiveX * dx; // current real value
                    prealBlock[x] = prealBlockStart[x] = startX;
               }

                //limit loop
                for (int lo = 0; lo < limit; lo++) {

                    int doneCount=0;
                    #pragma omp simd aligned(ppdata,prealBlock,pimagBlock,prealBlockStart) reduction(+:doneCount)
                    for (int x = 0; x < BLOCK; x++) {
                        float r2 = prealBlock[x] * prealBlock[x];
                        float i2 = pimagBlock[x] * pimagBlock[x];

                        if (r2 + i2 > 4.0f) {
                            doneCount++;
                            continue;
                        }

                        pimagBlock[x] = 2.0f * prealBlock[x] * pimagBlock[x] + startY;
                        prealBlock[x] = r2 - i2 + prealBlockStart[x];
                        //new effective X via pointer arithmetics
                        ppdata[x]++;

                    }
                    if(doneCount==BLOCK)
                        break;
                }
            }
        }
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