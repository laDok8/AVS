/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  LADISLAV DOKOUPIL <xdokou14@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    DATE
 **/

#include <iostream>
#include <math.h>
#include <limits>

#include "tree_mesh_builder.h"

TreeMeshBuilder::TreeMeshBuilder(unsigned gridEdgeSize)
    : BaseMeshBuilder(gridEdgeSize, "Octree")
{

}

unsigned TreeMeshBuilder::marchCubes(const ParametricScalarField &field)
{
    // Suggested approach to tackle this problem is to add new method to
    // this class. This method will call itself to process the children.
    // It is also strongly suggested to first implement Octree as sequential
    // code and only when that works add OpenMP tasks to achieve parallelism.

    return marchSegment(Vec3_t<float>(0,0,0), mGridSize, field);
}

unsigned TreeMeshBuilder::marchSegment(const Vec3_t<float> &cubeOffset,const unsigned int cubeLen, const ParametricScalarField &field)
{
    unsigned totalTriangles = 0;

    // 2. Loop over each coordinate in the 3D grid.
    for(size_t i = 0; i < 8; ++i)
    {
        const float edgeLength = static_cast<float>(cubeLen) / 2.f;

        Vec3_t<float> cubeStart( cubeOffset.x + ((i % 2) * edgeLength),
                                  cubeOffset.y + (((i / 2) % 2)*edgeLength),
                                  cubeOffset.z + ((i / (2*2))*edgeLength));

        const float add = edgeLength * mGridResolution / 2.f;
        Vec3_t<float> midPoint(cubeStart.x * mGridResolution  + add, cubeStart.y * mGridResolution + add, cubeStart.z * mGridResolution  + add);

        float actualISO = evaluateFieldAt(midPoint,field);
        if(actualISO > mIsoLevel + edgeLength*mGridResolution*sqrt(3)/2.0f) {
            continue;
        }


        if(cubeLen>32){
            totalTriangles += marchSegment(cubeStart,cubeLen/2,field);
        }
        else{
            totalTriangles += marchSmallCube(cubeStart, cubeLen/2, field);
        }


    }

    return totalTriangles;
}

//build small cubes
unsigned TreeMeshBuilder::marchSmallCube(const Vec3_t<float> &pos, const unsigned int cubeLen, const ParametricScalarField &field){
    size_t totalCubesCount = cubeLen*cubeLen*cubeLen;

    unsigned totalTriangles = 0;

    for(size_t i = 0; i < totalCubesCount; ++i)
    {
        // 3. Compute 3D position in the grid.
        Vec3_t<float> cubeOffset( pos.x + (i % cubeLen),
                                  pos.y + ((i / cubeLen) % cubeLen),
                                  pos.z + (i / (cubeLen*cubeLen)));

        // 4. Evaluate "Marching Cube" at given position in the grid and
        //    store the number of triangles generated.
        totalTriangles += buildCube(cubeOffset, field);
    }
    return totalTriangles;

}

float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    // 1. Store pointer to and number of 3D points in the field
    //    (to avoid "data()" and "size()" call in the loop).
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    // 2. Find minimum square distance from points "pos" to any point in the
    //    field.
    for(unsigned i = 0; i < count; ++i)
    {
        float distanceSquared  = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared       += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared       += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        // Comparing squares instead of real distance to avoid unnecessary
        // "sqrt"s in the loop.
        value = std::min(value, distanceSquared);
    }

    // 3. Finally, take square root of the minimal square distance to get the real distance
    return sqrt(value);
}

void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    mTriangles.push_back(triangle);
}
