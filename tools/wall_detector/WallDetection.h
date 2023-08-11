#ifndef WALLDETECTION_H
#define WALLDETECTION_H

bool isNormallyDistributed(const std::vector<double>& data, double significance_level = 0.05);

Eigen::Vector4d findMinimizingPlane(const vector<Eigen::Vector3d>& points);

double angleBetweenPlanes(Eigen::Vector3d normalVector);

bool wallDetector(vector <Eigen::Vector3d>& points);

#endif //WALLDETECTION_H