#include<iostream>
#include<fstream>
#include<iomanip>
#include"adjustment.h"

int main(){
    system("chcp 65001 > nul");
    std::vector<GroundPoint> points;
    readfile(points);

    loadAllExteriorOrientations();
    std::cout<<"已加载 "<<imageData.size()<<" 个影像的外方位元素"<<std::endl;

    double f = 105.2;
    int failCount = 0;

    for(size_t i = 0; i < points.size(); i++){
        if(!adj(points[i], f, i)){
            failCount++;
        }
    }

    std::ofstream outfile("result.txt");
    for(size_t i = 0; i < points.size(); i++){
        outfile << std::fixed << std::setprecision(4)
                << points[i].X << " " << points[i].Y << " " << points[i].Z << std::endl;
    }
    outfile.close();

    std::cout << "计算完成，结果已保存到 result.txt" << std::endl;
    std::cout << "共处理 " << points.size() << " 个点" << std::endl;
    std::cout << "平差失败点: " << failCount << " 个" << std::endl;
    return 0;
}
