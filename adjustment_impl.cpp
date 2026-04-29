#include<iostream>
#include<fstream>
#include<string>
#include<limits>
#include<cmath>
#include"adjustment.h"

std::map<int, ExteriorOrientation> imageData;

void readfile(std::vector<GroundPoint>& points){
    std::ifstream file("Data.pts");
    if (!file.is_open()){
        std::cout<<"无法打开Data.pts 文件！"<<std::endl;
        return;
    }
    std::string line;
    for(int i=0;i<6;i++){
        getline(file, line);
    }

    while(!file.eof()){
        GroundPoint pt;
        std::string pointID;
        int attrib;
        if(!(file>>pointID>>pt.X>>pt.Y>>pt.Z>>attrib)) break;
        pt.X = std::round(pt.X / 10.0) * 10.0;
        pt.Y = std::round(pt.Y / 10.0) * 10.0;
        pt.Z = std::round(pt.Z / 10.0) * 10.0;
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        int numImages;
        file>>numImages;
        for(int j=0; j<numImages; j++){
            int ImageID, groupID;
            double x, y;
            file>>ImageID>>x>>y>>groupID;
            ImagePoint P;
            P.imageID = ImageID;
            P.x = x;
            P.y = y;
            pt.imagePoints.push_back(P);
        }
        points.push_back(pt);
    }
    file.close();
}

void loadAllExteriorOrientations(){
    std::ifstream file("Data.pht");
    if (!file.is_open()){
        std::cout<<"Data.pht无法打开"<<std::endl;
        return;
    }
    std::string line;
    for(int i=0; i<3; i++){
        getline(file, line);
    }

    int totalImage;
    file>>totalImage;
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    for(int i=0; i<totalImage; i++){
        int id;
        file>>id;
        ExteriorOrientation eo;
        file>>eo.Xs>>eo.Ys>>eo.Zs>>eo.Phi>>eo.Omega>>eo.Kappa;
        imageData[id] = eo;
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    file.close();
}

void find_Image_ExteriorOrientation(ExteriorOrientation& eo, int ImageID){
    if(imageData.find(ImageID) != imageData.end()){
        eo = imageData[ImageID];
    }
}

Eigen::Matrix3d getRotationMatrix(const ExteriorOrientation& eo){
    Eigen::Matrix3d R;
    R(0, 0) = std::cos(eo.Phi)*std::cos(eo.Kappa)-std::sin(eo.Phi)*std::sin(eo.Omega)*std::sin(eo.Kappa);
    R(0, 1) = -std::cos(eo.Phi)*std::sin(eo.Kappa)-std::sin(eo.Phi)*std::sin(eo.Omega)*std::cos(eo.Kappa);
    R(0, 2) = -std::sin(eo.Phi)*std::cos(eo.Omega);
    R(1, 0) = std::cos(eo.Omega)*std::sin(eo.Kappa);
    R(1, 1) = std::cos(eo.Omega)*std::cos(eo.Kappa);
    R(1, 2) = -std::sin(eo.Omega);
    R(2, 0) = std::sin(eo.Phi)*std::cos(eo.Kappa)+std::cos(eo.Phi)*std::sin(eo.Omega)*std::sin(eo.Kappa);
    R(2, 1) = -std::sin(eo.Phi)*std::sin(eo.Kappa) + std::cos(eo.Phi)*std::sin(eo.Omega)*std::cos(eo.Kappa);
    R(2,2) = std::cos(eo.Phi)*std::cos(eo.Omega);
    return R;
}

bool adj(GroundPoint& gp, double f, int pointIndex){
    const int maxIter = 100;
    const double threshold = 1e-6;
    std::vector<bool> valid(gp.imagePoints.size(), true);
    bool hasOutlier = false;

    

    for(int iter = 0; iter < maxIter; iter++){
        std::vector<int> validIdx;
        for(size_t i = 0; i < gp.imagePoints.size(); i++){
            if(valid[i]) validIdx.push_back(i);
        }

        int n = validIdx.size() * 2;
        if(n < 6) return false;

        Eigen::MatrixXd B(n, 3);
        Eigen::VectorXd l(n);

        for(size_t k = 0; k < validIdx.size(); k++){
            size_t i = validIdx[k];
            ExteriorOrientation eo;
            find_Image_ExteriorOrientation(eo, gp.imagePoints[i].imageID);
            Eigen::Matrix3d R = getRotationMatrix(eo);

            double X = gp.X - eo.Xs;
            double Y = gp.Y - eo.Ys;
            double Z = gp.Z - eo.Zs;

            double a1 = R(0,0), a2 = R(0,1), a3 = R(0,2);
            double b1 = R(1,0), b2 = R(1,1), b3 = R(1,2);
            double c1 = R(2,0), c2 = R(2,1), c3 = R(2,2);

            double X_bar = a1*X + b1*Y + c1*Z;
            double Y_bar = a2*X + b2*Y + c2*Z;
            double Z_bar = a3*X + b3*Y + c3*Z;

            double denom = Z_bar * Z_bar;

            B(2*k, 0) = -f * (a1*Z_bar - a3*X_bar) / denom;
            B(2*k, 1) = -f * (b1*Z_bar - b3*X_bar) / denom;
            B(2*k, 2) = -f * (c1*Z_bar - c3*X_bar) / denom;

            B(2*k+1, 0) = -f * (a2*Z_bar - a3*Y_bar) / denom;
            B(2*k+1, 1) = -f * (b2*Z_bar - b3*Y_bar) / denom;
            B(2*k+1, 2) = -f * (c2*Z_bar - c3*Y_bar) / denom;

            double x_img = gp.imagePoints[i].x ;
            double y_img = gp.imagePoints[i].y ;

            double x_calc = -f * X_bar / Z_bar;
            double y_calc = -f * Y_bar / Z_bar;

            l(2*k) = x_img - x_calc;
            l(2*k+1) = y_img - y_calc;
        }

        Eigen::Vector3d delta = (B.transpose()*B).ldlt().solve(B.transpose()*l);

        gp.X += delta(0);
        gp.Y += delta(1);
        gp.Z += delta(2);

        Eigen::VectorXd V = B * delta - l;
        double sigma0 = (n > 3) ? std::sqrt((V.transpose()*V)(0) / (n - 3)) : 0;

        for(size_t k = 0; k < validIdx.size(); k++){
            double v = std::sqrt(V(2*k)*V(2*k) + V(2*k+1)*V(2*k+1));
            if(v > 3.0 * sigma0 && sigma0 > 0){
                valid[validIdx[k]] = false;
                if(!hasOutlier){
                    std::cout << "点 " << pointIndex << " 检测到粗差，剔除在影像 " << gp.imagePoints[validIdx[k]].imageID <<"上的观测"<< std::endl;
                    hasOutlier = true;
                } else {
                    std::cout << "  继续剔除影像 " << gp.imagePoints[validIdx[k]].imageID << std::endl;
                }
            }
        }

        if(delta.norm() < threshold) break;
    }

    return true;
}
