//
// Created by Wei Chen on 22-6-29.
//

#include "MarchingCubes.hpp"
#include "Utils.h"

#include <spdlog/spdlog.h>
#include <igl/write_triangle_mesh.h>
#include <boost/progress.hpp>

namespace MC {

MarchingCubes::MarchingCubes(int p_nX, int p_nY, int p_nZ, unsigned short p_iso_value, const std::string &data_file) :
        nX(p_nX), nY(p_nY), nZ(p_nZ), iso_value(p_iso_value)
{
    data.resize(nZ);
    for (int zI=0; zI<nZ; ++zI) {
        data[zI].resize(nY);
        for (int yI=0; yI<nY; ++yI) {
            data[zI][yI].resize(nX);
        }
    }

    if(!read_data(data_file, data)) {
        spdlog::error("error on reading data!");
        exit(-1);
    }

    nV.resize(nX-1);
    nF.resize(nX-1);
    triV.resize(nX-1);
    triF.resize(nX-1);
    for (int xI=0; xI < nX-1; ++xI) {
        nV[xI] = 0;
        nF[xI] = 0;
        triV[xI].resize(nV[xI], 3);
        triF[xI].resize(nF[xI], 3);
    }

    spdlog::info("MarchingCubes constructed");
}

void MarchingCubes::solve() {
    boost::progress_display show_progress(nX);

#pragma omp parallel for
    for (int i=0; i<nX-1; i++) {
        for (int j=0; j<nY-1; j++) {
            for (int k=0; k<nZ-1; k++) {
                Eigen::Matrix<double, 8, 3> V;
                V <<   i,   j,   k,
                     i+1,   j,   k,
                     i+1, j+1,   k,
                       i, j+1,   k,
                       i,   j, k+1,
                     i+1,   j, k+1,
                     i+1, j+1, k+1,
                       i, j+1, k+1;
                Eigen::Vector<unsigned short, 8> val;
                val << data[  k][  j][  i],
                       data[  k][  j][i+1],
                       data[  k][j+1][i+1],
                       data[  k][j+1][  i],
                       data[k+1][  j][  i],
                       data[k+1][  j][i+1],
                       data[k+1][j+1][i+1],
                       data[k+1][j+1][  i];
                one_cell(i, V, val);
            }
        }

        ++show_progress;
    }

}

void MarchingCubes::one_cell(int xI, const Eigen::Matrix<double, 8, 3> &V, const Eigen::Vector<unsigned short, 8> &val){
    int cubeIndex = 0;
    Eigen::Matrix<double, 12, 3> vertlist;
    // compute index
    if (val[0] < iso_value) cubeIndex |= 1;
    if (val[1] < iso_value) cubeIndex |= 2;
    if (val[2] < iso_value) cubeIndex |= 4;
    if (val[3] < iso_value) cubeIndex |= 8;
    if (val[4] < iso_value) cubeIndex |= 16;
    if (val[5] < iso_value) cubeIndex |= 32;
    if (val[6] < iso_value) cubeIndex |= 64;
    if (val[7] < iso_value) cubeIndex |= 128;

    // no intersected points
    if (edgeTable[cubeIndex] == 0 || edgeTable[cubeIndex] == 255) {
        return;
    }
    // compute intersected points using interpolation
    if (edgeTable[cubeIndex] & 1) {
        vertlist.row(0) = VertexInterp(iso_value, V.row(0), V.row(1), val(0), val(1));
    }
    if (edgeTable[cubeIndex] & 2) {
        vertlist.row(1) = VertexInterp(iso_value, V.row(1), V.row(2), val(1), val(2));
    }
    if (edgeTable[cubeIndex] & 4) {
        vertlist.row(2) = VertexInterp(iso_value, V.row(2), V.row(3), val(2), val(3));
    }
    if (edgeTable[cubeIndex] & 8) {
        vertlist.row(3) = VertexInterp(iso_value, V.row(3), V.row(0), val(3), val(0));
    }
    if (edgeTable[cubeIndex] & 16) {
        vertlist.row(4) = VertexInterp(iso_value, V.row(4), V.row(5), val(4), val(5));
    }
    if (edgeTable[cubeIndex] & 32) {
        vertlist.row(5) = VertexInterp(iso_value, V.row(5), V.row(6), val(5), val(6));
    }
    if (edgeTable[cubeIndex] & 64) {
        vertlist.row(6) = VertexInterp(iso_value, V.row(6), V.row(7), val(6), val(7));
    }
    if (edgeTable[cubeIndex] & 128) {
        vertlist.row(7) = VertexInterp(iso_value, V.row(7), V.row(4), val(7), val(4));
    }
    if (edgeTable[cubeIndex] & 256) {
        vertlist.row(8) = VertexInterp(iso_value, V.row(0), V.row(4), val(0), val(4));
    }
    if (edgeTable[cubeIndex] & 512) {
        vertlist.row(9) = VertexInterp(iso_value, V.row(1), V.row(5), val(1), val(5));
    }
    if (edgeTable[cubeIndex] & 1024) {
        vertlist.row(10) = VertexInterp(iso_value, V.row(2), V.row(6), val(2), val(6));
    }
    if (edgeTable[cubeIndex] & 2048) {
        vertlist.row(11) = VertexInterp(iso_value, V.row(3), V.row(7), val(3), val(7));
    }

    for (int i=0; triTable[cubeIndex][i] != -1; i+=3) {
        triV[xI].conservativeResize(nV[xI]+3, 3);
        triF[xI].conservativeResize(nF[xI]+1, 3);
        triV[xI].row(nV[xI]) = vertlist.row(triTable[cubeIndex][i  ]);
        triV[xI].row(nV[xI]+1) = vertlist.row(triTable[cubeIndex][i+1]);
        triV[xI].row(nV[xI]+2) = vertlist.row(triTable[cubeIndex][i+2]);
        triF[xI].row(nF[xI]) << nV[xI], nV[xI]+1, nV[xI]+2;
        nV[xI] += 3;
        nF[xI] += 1;
    }

}

void MarchingCubes::save(const std::string &savePath) {
    spdlog::info("saving result to obj file");

    int nV_ = 0;
    int nF_ = 0;
    for (int xI=0; xI<nX-1; ++xI) {
        nV_ += nV[xI];
        nF_ += nF[xI];
    }
    Eigen::MatrixXd triV_(nV_, 3);
    Eigen::MatrixXi triF_(nF_, 3);
    int tmpV = 0;
    int tmpF = 0;
    for (int xI=0; xI<nX-1; ++xI) {
        triV_(Eigen::seq(tmpV, tmpV+nV[xI]-1), Eigen::all) = triV[xI];
        triF_(Eigen::seq(tmpF, tmpF+nF[xI]-1), Eigen::all) = triF[xI].array() + tmpV;
        tmpV += nV[xI];
        tmpF += nF[xI];
    }
    igl::write_triangle_mesh(savePath, triV_, triF_);
}

} // MC
