
#ifndef __DE_H_
#define __DE_H_
#include "func.h"
#include <time.h>

class RandomNumberGenerator {
public:
	RandomNumberGenerator() { srand(time(0)); }
	// get a random number in [low, high]
	double randomFBetween(double low, double high) {
		double d = (double)rand() / (double)RAND_MAX;
		return d * (high - low) + low;
	}
	int randomIBetween(int low, int high) {
		int v = rand() % (high - low + 1);
		return v + low;
	}
};

class DE
{
    public:
        DE(int p, int dp, int di, int gmax, vector<string> trainFilePathSet, BaseFunc* func, vector<double>& initParams, vector<double>& paramsShockRatio, vector<string> paramsNameVec);

        void init(string filePath);
        vector<double> startDE();
        void saveBestParams();
        void saveBestResult();

    public:
        int p;      // population size
        int dp;      // parameters num(dimension of parameters)
        int di;     // dimension of input
        int gmax;   // the number of maximum generations
        BaseFunc* func;     // fitness function
        vector<vector<double>> paramsPopulation;  // the set of population of p individuals
        vector<double> error;           // the error of fitness func, the dimension is p
        vector<double> bestParams;      // best params, dimension is d
        double bestError;               // the error of best params

        vector<double> initParams;      // initial parameters
        vector<double> paramsShockRatio;        // the maximum ratio of shock parameters
        double initError;               // initial error

        vector<double> bestResult;
    public:
        vector<vector<vector<double>>> trainDatasetVec;    // train dataset, dimension is d + 1,
                                        // first p value is input, last value is real observations
        vector<double> trainDataSetMaxObserve;
        vector<string> trainDatasetPathVec;

        vector<string> paramsNameVec;
    private:
        void readTrainDataset(vector<string> filePathSet);
        void checkParams(double& param, int dimension);
        double calcuError(vector<double>& params);
        double computeBest();

        FILE *optFileStream = NULL;
        string optFilePath;
    public:
        void setOptFile(string filePath);
        void flushFile();
        void writeFile(int generationNum, double error);
        void closeFile();
};

#endif // __DE_H_
