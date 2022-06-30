//
// Created by Wei Chen on 22-6-29.
//

#ifndef MARCHINGCUBES_UTILS_H
#define MARCHINGCUBES_UTILS_H

#include <Eigen/Eigen>
#include <spdlog/spdlog.h>

namespace MC{

bool read_data(const std::string &filePath, std::vector<std::vector<std::vector<unsigned short>>> &data){
    FILE* in = fopen(filePath.c_str(), "rb");
    if (!in) {
        return false;
    }

    int Min = 0x3f3f3f3f, Max = -0x3f3f3f3f;
    for (int zI=0; zI<data.size(); ++zI) {
        for (int yI=0; yI<data[zI].size(); ++yI) {
            fread(data[zI][yI].data(), sizeof(data[zI][yI][0]), data[zI][yI].size(), in);
            for (const auto item : data[zI][yI]) {
                Min = std::min(Min, (int)item);
                Max = std::max(Max, (int)item);
            }
        }
    }
    spdlog::info("Min value in data: {}", Min);
    spdlog::info("Max value in data: {}", Max);

    fclose(in);
    return true;
}

Eigen::RowVector3d VertexInterp(unsigned short iso_value, const Eigen::RowVector3d &p1, const Eigen::RowVector3d &p2,
                                unsigned short val1, unsigned short val2) {
    assert(iso_value >= std::min(val1, val2) && iso_value <= std::max(val1, val2));
    if (val1 == val2) {
        return p1;
    }
    double scalar = 1.0 * (iso_value - val1) / (val2 - val1);
    return p1 + scalar * (p2 - p1);
}

} // namespace MC

#endif //MARCHINGCUBES_UTILS_H
