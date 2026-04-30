#ifndef ADJUSTMENT_H
#define ADJUSTMENT_H

#include<vector>
#include<map>
#include"Eigen/Dense"

struct ImagePoint{
    int imageID;
    double x, y;
};

struct GroundPoint{
    double X, Y, Z;
    std::vector<ImagePoint> imagePoints;
};

struct ExteriorOrientation{
    double Xs, Ys, Zs;
    double Phi, Omega, Kappa;
};

extern std::map<int, ExteriorOrientation> imageData;

void readfile(std::vector<GroundPoint>& points);
void loadAllExteriorOrientations();
void find_Image_ExteriorOrientation(ExteriorOrientation& eo, int ImageID);
Eigen::Matrix3d getRotationMatrix(const ExteriorOrientation& eo);
bool adj(GroundPoint& gp, double f, int pointIndex);

#endif
