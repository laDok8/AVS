/**
 * @file BatchMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
{
    data = (int *)(malloc(height * width * sizeof(int)));
    realBlock = (float *)(malloc(BLOCK * sizeof(float)));
    imagBlock = (float *)(malloc(BLOCK * sizeof(float)));
    realBlockStart = (float *)(malloc(BLOCK * sizeof(float)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
    free(data);
    data = NULL;
    free(realBlock);
    realBlock = NULL;
    free(imagBlock);
    imagBlock = NULL;
    free(realBlockStart);
    realBlockStart = NULL;
}


int * BatchMandelCalculator::calculateMandelbrot () {
    int *pdata = data;

    //TODO: staci height/2
    //y block
    for (int ty = 0; ty < height / BLOCK; ty++) {
        //y block
        for (int tx = 0; tx < width / BLOCK; tx++) {

            for (int y = 0; y < BLOCK; y++) {
                const int effectiveY = ty * BLOCK + y;
                if(effectiveY > height) {
                    break;
                }

                //init values
                const float startY = y_start + effectiveY * dy; // current imaginary value
                for (int x = 0; x < BLOCK; x++) {
                    imagBlock[x] = startY;
                    const int effectiveX = tx * BLOCK + x;
                    const float startX = x_start + effectiveX * dx; // current real value
                    realBlock[x] = realBlockStart[x] = startX;
               }

                //limit loop
                for (int lo = 0; lo < limit; lo++) {

                    for (int x = 0; x < BLOCK; x++) {
                        const int effectiveX = tx * BLOCK + x;
                        if(effectiveX > width) {
                            break;
                        }

                        float r2 = realBlock[x] * realBlock[x];
                        float i2 = imagBlock[x] * imagBlock[x];

                        if (r2 + i2 > 4.0f) {
                            continue;
                        }

                        imagBlock[x] = 2.0f * realBlock[x] * imagBlock[x] + startY;
                        realBlock[x] = r2 - i2 + realBlockStart[x];
                        pdata[effectiveY * width + effectiveX]++;

                    }
                }
            }
        }
    }
    return data;
}