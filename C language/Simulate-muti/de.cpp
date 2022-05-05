
#include <fstream>
#include <string>
#include <sstream>
#include <limits.h>
#include <float.h>
// #include <math.h>
#include <cmath>
#include <thread>
#include <omp.h>

#include "func.h"
#include "de.h"

#ifdef _WIN32
    #include <direct.h>
    #define SIZE_PT "llu"
#elif defined(__linux__)
    #include <sys/stat.h>
    #include <unistd.h>
    #define SIZE_PT "lu"
#endif

void makeDir(string dirPath)
{
    #ifdef _WIN32
        mkdir(dirPath.c_str());
    #elif defined(__linux__)
        mkdir(dirPath.c_str(), 0777);
    #endif
}

bool fileExist(string filePath)
{
    #ifdef _WIN32
        return access(filePath.c_str(), 0) == 0;
    #elif defined(__linux__)
        return access(filePath.c_str(), 0) == F_OK;
    #endif
    return false;
}

RandomNumberGenerator rng;

DE::DE(int p, int dp, int di, int gmax,
    vector<string> trainFilePathSet,
    BaseFunc* func,
    vector<double>& initParams, vector<double>& paramsShockRatio,
    vector<string> paramsNameVec
    )
    :
    p(p), dp(dp), di(di), gmax(gmax),
    paramsPopulation(p), error(p),
    func(func),
    trainDatasetVec(trainFilePathSet.size()), trainDataSetMaxObserve(trainFilePathSet.size()), trainDatasetPathVec(trainFilePathSet),
    initParams(initParams), paramsShockRatio(paramsShockRatio),
    paramsNameVec(paramsNameVec)
{
    // read train dataset
    this->readTrainDataset(trainFilePathSet);

    // generate initial populition
    paramsPopulation[0] = initParams;
    for(size_t i=1; i<p; i++){
        this->paramsPopulation[i].resize(dp);
        for(size_t j=0; j<dp; j++){
            this->paramsPopulation[i][j] = initParams[j]
                + rng.randomFBetween(-paramsShockRatio[j], paramsShockRatio[j])
                    * initParams[j];
        }
    }
    // calculate initial fitness error
    for(size_t i=0; i<p; i++){
        this->error[i] = calcuError(this->paramsPopulation[i]);
    }
    // init best params
    double initError = this->computeBest();
    this->initError = initError;
    printf("generation 0: error is: %g\n", this->initError);
}

double DE::calcuError (vector<double>& params)
{
    double totalError = 0;
    for(size_t j=0; j<this->trainDatasetVec.size(); j++){
        vector<vector<double>>& trainSet = this->trainDatasetVec[j];
        double maxObserv = this->trainDataSetMaxObserve[j];

        double error = 0;
        double sumRatio = 1.2;
        double logRatio = 0.05;

        for(size_t i=0; i<trainSet.size(); i++){
            vector<double> input;
            input.assign(
                trainSet[i].begin(),
                trainSet[i].end() - 1
            );
            double output = func->func(params, input, this->di);
    //        if(output == 0) continue;

         
            error += pow((trainSet[i].back()-output)/maxObserv, 2.0); //线性 RMS ERROR
         //   error +=  pow( (log10(trainSet[i].back()) - log10(output) )/ log10(maxObserv), 2.0); //对数RMS ERROR

         //   error = error / 2;


     /*       if(trainSet[i].back() < 2*10e-5){
                error += logRatio * pow( (log10(trainSet[i].back()) - log10(output) )/ log10(maxObserv), 2.0); //对数RMS ERROR
            }else{
                error += sumRatio * pow((trainSet[i].back()-output)/maxObserv, 2.0); //线性 RMS ERROR
            }
*/
        }
        error /= (double)trainSet.size();

        error = sqrt(error) * 100;


      //给每幅图加权

        // printf("%d: %.2f\n", j, error);
/*
 if(  j == 0  ){
            totalError += 1 * error;}
 else if(  j == 1  ){
            totalError += 1 * error;}
 else{
            totalError += error;}
*/
       totalError += error;
    }

    return (this->trainDatasetVec.size()==0 ? 0 : totalError / this->trainDatasetVec.size());

//   return totalError / 3;
}

double DE::computeBest()
{
    int idx = 0;
    double minError = this->error[0];
    for(size_t i=1; i<this->error.size(); i++){
        if(this->error[i] < minError){
            minError = this->error[i];
            idx = i;
        }
    }
    this->bestParams = this->paramsPopulation[idx];
    this->bestError = minError;
    return this->bestError;
}

void DE::readTrainDataset (vector<string> filePathSet)
{
    for(size_t i=0; i<filePathSet.size(); i++){
        string filePath = filePathSet[i];
        ifstream fin(filePath);
        if(! fin){
            printf("can not read train dataset at: %s\n", filePath.c_str());
            exit(0);
        }
        string buffer;
        double maxObserve = DBL_MIN;
        // eat first line
        getline(fin, buffer);

        stringstream ss;
        while(getline(fin, buffer)){
            ss.clear();
            ss << buffer;
            double tmp;
            vector<double> valueVec;
            while(ss >> tmp){
                valueVec.push_back(tmp);
            }

            if(! valueVec.empty() && valueVec.size() == this->di + 1){

                if(valueVec.back() <= 1E-11){
                    continue;
                } // 如果观测值小于1e-11，则忽略该组数据

                this->trainDatasetVec[i].push_back(valueVec);
                maxObserve = max(maxObserve, valueVec.back());
            }else{
                // 数据集格式有误

            }
        }
        this->trainDataSetMaxObserve[i] = maxObserve;

        printf("train set size of file %s: %" SIZE_PT ", max value: %g\n", filePath.c_str(), this->trainDatasetVec[i].size(), maxObserve);
        fin.close();
    }

}

void DE::checkParams (double& param, int dimension)
{
    if(dimension<0 || dimension>this->initParams.size()) return ;
    double maxParam = this->initParams[dimension] * (1 + this->paramsShockRatio[dimension]);
    double minParam = this->initParams[dimension] * (1 - this->paramsShockRatio[dimension]);
    if(param > maxParam){
        param = maxParam;
    }
    if(param < minParam){
        param = minParam;
    }
}
int getCpuNum()
{
    /*
	char buf[16] = {0};
	int num = 1;
	FILE *fp = popen("cat /proc/cpuinfo | grep processor | wc -l", "r");
	if (fp)
	{
		fread(buf, 1, sizeof(buf) - 1, fp);
		pclose(fp);
	}
	num = atoi(buf);
	num = (num <= 0 ? 1 : num);
	return num;
	*/
	return thread::hardware_concurrency();
}
vector<double> DE::startDE()
{
    // 迭代次数
    int threadsNum = ceil(getCpuNum() * 0.8);
    omp_set_num_threads(threadsNum);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        for(int loop=tid; loop<this->gmax; loop+=threadsNum){
            // 对每个个体
            for (int i=tid, cnt=0; cnt<p; i=(i+1)%p, cnt++ ) {

                vector<double> xr1 = this->paramsPopulation[rng.randomIBetween(0, this->p - 1)];
                vector<double> xr2 = this->paramsPopulation[rng.randomIBetween(0, this->p - 1)];
                vector<double> xr3 = this->paramsPopulation[rng.randomIBetween(0, this->p - 1)];

                vector<double> params_i = this->paramsPopulation[i];

                // mutation
                vector<double> v(this->dp); // 产生新个体1
                for (int j = 0; j < this->dp; j++) {   // 针对所有参数
                    double F = rng.randomFBetween(0.0, 1.0);

    //给每个参数附上不同的F

      /*   	    if( j == 3){
                   F = rng.randomFBetween(0.0, 5.0);
                            }
            else if( j == 13 || j == 14){
                        F = rng.randomFBetween(0.0, 10.0);
                    }
*/
                    v[j] = xr1[j] + F * (xr2[j] - xr3[j]);

                    if((j == 0 ||j == 1 ||j == 2 || j == 3 || j == 4 || j == 6 ||j == 7 || j == 9|| j == 8)
                        && v[j] < 0){
                        v[j] = params_i[j];}

                    if((j == 10)
                        && (v[j] < 0.4 || v[j] >0.465)){
                        v[j] = params_i[j];}
/*                     if((j >= 11 || j<=15 )
                        && (v[j] < 0.7 || v[j] >1.3)){
                        v[j] = params_i[j];} */
 /*
                    if((j == 5 ) && (v[j]<0.000122 ||v[j]>0.000128))  //alpha
                        {
                        v[j] = params_i[j];}

                    if((j == 6 ) && (v[j]<5.5e-21||v[j]>5.9e-21))  //beeta
                        {
                        v[j] = params_i[j];}

                    if((j == 7 ) && (v[j]<1e40||v[j]>1e50))   //gammma
                        {
                        v[j] = params_i[j];}

                    if((j == 9 ) && (v[j]<1.45e-12 ||v[j]>1.57e-12))   //qo
                        {
                        v[j] = params_i[j];}

                    if((j == 12 ) && (v[j]<0.85 ||v[j]>1.2))   //Eta
                    {
                        v[j] = params_i[j];}

                    if((j == 10 ) && (v[j]<0 ||v[j]>1E-6))  //dl
                        {
                        v[j] = params_i[j];}

                if((j == 11 ) && (v[j]<0 ||v[j]>1E-6))   //dw
                        {
                        v[j] = params_i[j];}

  */

            //		checkParams(v[j], j); //限制参数在初始值的一定范围内波动
                }
                // crossover
                vector<double> u(this->dp); // 产生新个体2
                for (int j = 0; j < this->dp; j++) {
                    double CR = rng.randomFBetween(0.0, 1.0);
                    double CRx = rng.randomFBetween(0.0, 1.0);
                    if (CRx < CR) {
                        u[j] = v[j]; // 变异后的参数
                    } else {
                        u[j] = params_i[j];  // 保留原来的参数
                    }
                }

    /*
    // 进化策略算法
                int pos = rng.randomIBetween(0, this->dp);
                vector<double> u(this->dp);
                u.assign(
                    this->paramsPopulation[i].begin()
                    , this->paramsPopulation[i].end()
                );
                u[pos] = u[pos] * (1 + 0.01 * rng.randomFBetween(-1.0, 1.0));
    */

                double fit_error = this->calcuError(u);
                if (fit_error < this->error[i]) {
                    #pragma omp critical
                    {
                        this->paramsPopulation[i] = u;
                        this->error[i] = fit_error;
                    }
                }
            }
            #pragma omp critical
            {
                this->computeBest();
                // 写入文件
                this->writeFile(loop, this->bestError);
            }
            printf("thread#%d: generation %d: error is %g\n", tid, loop, this->bestError);


            if(loop == 5000
                || loop == 10000
                || loop == 25000
                || loop == 35000
                || loop == 50000
            )
            {

                this->saveBestParams();
                this->saveBestResult();
            }
        }
    }
	return this->bestParams;
}


void DE::saveBestParams()
{
    // string optFilePath = "params.txt";
    char fileName[200];
    snprintf(fileName, 200, "./params/params-%g.txt", this->bestError);
    string optFilePath = fileName;

    string dirPath = optFilePath.substr(0, optFilePath.find_last_of("/"));
    printf("dir: %s\n", dirPath.c_str());
    /*
    if(0 != access(dirPath.c_str(), 0)){
        mkdir(dirPath.c_str());
    }*/
    if(! fileExist(dirPath)){
        makeDir(dirPath);
    }

    ofstream fout(optFilePath, ios::app);
    if(! fout){
        printf("无法保存文件: %s\n", optFilePath.c_str());
        exit(2);
    }
    fout << "init params\n" << "error: " << this->initError << endl;
    for(size_t i=0; i<this->initParams.size(); i++){
       // fout << this->paramsNameVec[i] << ": " << this->initParams[i] << endl; //输出参数前面带para名字
        fout << this->initParams[i] << endl;//输出参数前面带para名字
    }
    fout << endl;
    fout << "new params\n" << "error: " << this->bestError << endl;
    for(size_t i=0; i<this->bestParams.size(); i++){
        // fout << this->paramsNameVec[i] << ": " << this->bestParams[i] << endl;
        fout << this->bestParams[i] << endl;
    }
    fout << "****************************************\n\n";
    fout.close();

    printf("params save to file: %s\n", optFilePath.c_str());
}

void DE::saveBestResult()
{
    vector<vector<double>> bestResultVec;
    bestResultVec.resize(this->trainDatasetVec.size());

    for(size_t j=0; j<this->trainDatasetVec.size(); j++){
        vector<vector<double>>& trainSet = this->trainDatasetVec[j];

        for(size_t i=0; i<trainSet.size(); i++){
            vector<double> input;
            input.assign(
                trainSet[i].begin(),
                trainSet[i].end() - 1
            );
            double output = func->func(this->bestParams, input, this->di);
            bestResultVec[j].push_back(output);
        }
    }
    // printf("best result size: %lld\n", bestResultVec.size());

    for(size_t k=0; k<bestResultVec.size(); k++){
        vector<double>& bestResult = bestResultVec[k];

        string dataFilePath = this->trainDatasetPathVec[k];
        string fileName = dataFilePath.substr(0, dataFilePath.find_last_of("."));
        fileName = fileName.substr(fileName.find_last_of("/")+1);

        char tmpChar[200];
        snprintf(tmpChar, 200, "./output/result-%g-%s.txt", this->bestError, fileName.c_str());
        printf("result file: %s\n", tmpChar);

        string optFilePath = string(tmpChar);
        string dirPath = optFilePath.substr(0, optFilePath.find_last_of("/"));
        printf("dir: %s\n", dirPath.c_str());
        /*if(0 != access(dirPath.c_str(), 0)){
            mkdir(dirPath.c_str());
        }*/
        if(! fileExist(dirPath)){
            makeDir(dirPath);
        }

        ofstream fout(tmpChar);
        if(! fout){
            printf("不能保存结果文件: %s\n", tmpChar);
            exit(2);
        }

        for(size_t i=0; i<bestResult.size(); i++){
            for(size_t j=0; j<this->trainDatasetVec[k][i].size(); j++){
                fout << this->trainDatasetVec[k][i][j] << " ";
            }
            fout << bestResult[i] << endl;
        }
        // fout << "****************************************\n\n";
        fout.close();
    }
}


void DE::setOptFile(string filePath)
{
    this->optFilePath = filePath;
    this->optFileStream = fopen(filePath.c_str(), "a+");
    if(this->optFileStream == NULL){
        printf("无法写入文件%s\n", filePath.c_str());
        exit(0);
    }
}

void DE::flushFile()
{
    /*
    if(this->optFileStream != NULL){
        fclose(this->optFileStream);
    }
    this->optFileStream = fopen(this->optFilePath.c_str(), "a+");
    if(this->optFileStream == NULL){
        printf("无法写入文件%s\n", this->optFilePath.c_str());
    }
    */
    fflush(this->optFileStream);
}

void DE::writeFile(int generationNum, double error)
{
    fprintf(this->optFileStream, "%d %g\n", generationNum, error);

    // 每隔100代刷新一次文件
    if(generationNum % 5 == 0){
        this->flushFile();
    }
}

void DE::closeFile()
{
    if(this->optFileStream != NULL){
        fclose(this->optFileStream);
    }
}
