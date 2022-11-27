/**
 * @file    tree_mesh_builder.cpp
 *
 * @author  LADISLAV DOKOUPIL <xdokou14@stud.fit.vutbr.cz>
 *
 * @brief   Parallel Marching Cubes implementation using OpenMP tasks + octree early elimination
 *
 * @date    DATE
 **/

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
    unsigned triangles = 0;
    #pragma omp parallel default(none) shared(field, triangles)
    #pragma omp single nowait
        triangles = marchSegment(Vec3_t<float>(0, 0, 0), mGridSize, field);
    return triangles;
}

unsigned TreeMeshBuilder::marchSegment(const Vec3_t<float> &cubeOffset, const unsigned cubeLen,
                                       const ParametricScalarField &field)
{
    //empty
    const float add = cubeLen * mGridResolution / 2.f;
    Vec3_t<float> midPoint(cubeOffset.x * mGridResolution + add, cubeOffset.y * mGridResolution + add,
                           cubeOffset.z * mGridResolution + add);
    float actualISO = evaluateFieldAt(midPoint, field);
    if(actualISO > mIsoLevel + cubeLen * mGridResolution * sqrt(3) / 2.0f){
        return 0;
    }

    //cutoff
    if(cubeLen < 4){
        return marchSmallCube(cubeOffset, cubeLen, field);
    }


    unsigned totalTriangles = 0;
    for (size_t i = 0; i < 8; ++i) {
        #pragma omp task default(none) shared(cubeOffset, cubeLen, field, totalTriangles) firstprivate(i)
        {
            Vec3_t<float> cubeStart(cubeOffset.x + ((i % 2) * (cubeLen / 2)),
                                    cubeOffset.y + (((i / 2) % 2) * (cubeLen / 2)),
                                    cubeOffset.z + ((i / (2 * 2)) * (cubeLen / 2)));

            const unsigned triangles = marchSegment(cubeStart, cubeLen / 2, field);
            #pragma omp atomic update
            totalTriangles += triangles;
        }
    }

    #pragma omp taskwait
    return totalTriangles;
}

//build small cubes
unsigned TreeMeshBuilder::marchSmallCube(const Vec3_t<float> &pos, const unsigned cubeLen,
                                         const ParametricScalarField &field)
{
    size_t totalCubesCount = cubeLen * cubeLen * cubeLen;
    unsigned totalTriangles = 0;

    for (size_t i = 0; i < totalCubesCount; ++i) {
        Vec3_t<float> cubeOffset(pos.x + (i % cubeLen),
                                 pos.y + ((i / cubeLen) % cubeLen),
                                 pos.z + (i / (cubeLen * cubeLen)));

        totalTriangles += buildCube(cubeOffset, field);
    }
    return totalTriangles;

}

    float TreeMeshBuilder::evaluateFieldAt(const Vec3_t<float> &pos, const ParametricScalarField &field)
{
    const Vec3_t<float> *pPoints = field.getPoints().data();
    const unsigned count = unsigned(field.getPoints().size());

    float value = std::numeric_limits<float>::max();

    for (unsigned i = 0; i < count; ++i) {
        float distanceSquared = (pos.x - pPoints[i].x) * (pos.x - pPoints[i].x);
        distanceSquared += (pos.y - pPoints[i].y) * (pos.y - pPoints[i].y);
        distanceSquared += (pos.z - pPoints[i].z) * (pos.z - pPoints[i].z);

        value = std::min(value, distanceSquared);
    }

    return sqrt(value);
}

    void TreeMeshBuilder::emitTriangle(const BaseMeshBuilder::Triangle_t &triangle)
{
    #pragma omp critical
    mTriangles.push_back(triangle);
}
