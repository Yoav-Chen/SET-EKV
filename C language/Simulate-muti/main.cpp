//..........................................................................
// Differential evolution algorithm for model parameters optimization
//
// Adapt to windows and Linux systems
//
// Written by Mr. Yuanke Zhang and Mr. Wenjie Wang
// Copyright (c) belongs to Mr. Yuanke Zhang and Mr. Wenjie Wang
//
//  Last edit: 2021.9.22
//..........................................................................
#include "de.h"
#include "func.h"
#include "set_model.h"
#include  <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>


#include "de.cpp"
#include "set_model.cpp"
#include "set_final.cpp"

using namespace std;

pair< vector<double>, vector<double> > readParams(string filePath, vector<string>& paramsNameVec)
{
    vector<double> params;
    vector<double> paramsShockRatio;

    ifstream fin(filePath);
    if(! fin){
        printf("can not read params from: %s\n", filePath.c_str());
        exit(0);
    }
    string buffer;
    int line = 0;
    stringstream ss;
    while (getline(fin, buffer))
    {
        line ++;
        ss.clear();

        // 获取双引号“”内的参数名
        string nameBuf = "";
        int startT = buffer.find_first_of("\"");
        int endT = buffer.find_last_of("\"");

        // 兼容无双引号、只有一个双引号的情况
        if(startT != string::npos && endT != string::npos && startT!=endT){
            nameBuf = buffer.substr(startT+1, endT-startT-1);
        }
        if(startT == endT && startT != string::npos){
            printf("输入的参数文件: 《%s》第%d行格式有误，未设置正确的双引号\n", filePath.c_str(), line);
            exit(0);
        }
        // 没有参数名时，设置默认的参数名
        if(nameBuf.empty()){
            nameBuf = "param" + to_string(line);
        }

        // 设置buffer为剩余的参数值
        buffer = buffer.substr(endT+1);

        ss << buffer;
        // printf("%s\n", buffer.c_str());

        vector<double> tmpVec;
        double tmp;
        while(ss>>tmp){
            tmpVec.push_back(tmp);
            // printf("%g\n", tmp);
        }
        if(tmpVec.size() != 2){
            printf("输入的参数文件:%s第%d行格式有误\n", filePath.c_str(), line);
            exit(0);
        }
        params.push_back(tmpVec[0]);
        paramsShockRatio.push_back(tmpVec[1]);
        paramsNameVec.push_back(nameBuf);

        // printf("%g\n", tmp);
    }

    fin.close();

    return make_pair(params, paramsShockRatio);
}

int main()
{

    vector<string> paramsNameVec;

    pair< vector<double>, vector<double> > paramsPair = readParams("./input/init-paramss-set_final.txt", paramsNameVec);
    MatlabFunc matlabFunc;
    // string filePath = "./dataset/dataset.txt";
    vector<string> filePathSet =
    {
    //    "./dataset/dataset_1010_153K_VB0V.txt",
    //    "./dataset/dataset_1010_153K_Vd50mV.txt",
    //    "./dataset/dataset_1010_153K_Vd1.8V.txt"
        // "./dataset/dataset_1010_77K_Vd50mV_step5mV.txt",
        // "./dataset/dataset_1010_77K_Vd1.8V_step5mV.txt",
        // "./dataset/dataset_1010_77K_Vb0V_step50mV.txt"
        "./dataset/Id-Vg-Vd=0.002V.txt",
        "./dataset/Id-Vg-Vd=0.006V.txt"
    };

    DE de(
        2000                             // 种群大小
        , (int)paramsPair.first.size()          // 参数个数
        , 2                                    // 输入个数 (VG,VD,W,L,etc.)
        , 500                       // 最大迭代次数
        , filePathSet                              // 数据集文件路径
        , &matlabFunc                           // 待拟合的函数
        , paramsPair.first                      // 初始参数值
        , paramsPair.second                     // 参数波动范围
        , paramsNameVec                         // 参数名
    );

    // 设置写入文件
    de.setOptFile("generation-error.txt");

    de.startDE();
    de.saveBestParams();
    de.saveBestResult();

        // 关闭文件
    de.closeFile();

    return 0;
}

