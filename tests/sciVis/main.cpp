//
// Created by Wei Chen on 22-6-29.
//

#include "MarchingCubes.hpp"

int main() {
    MC::MarchingCubes mc = MC::MarchingCubes(512, 512, 507, 1000,
                                             "/home/cw/MyCode/MarchingCubes/input/visualization-data/raw_data/cbct_sample_z=507_y=512_x=512.raw");

    mc.solve();
    mc.save("/home/cw/MyCode/MarchingCubes/output/result.obj");
}