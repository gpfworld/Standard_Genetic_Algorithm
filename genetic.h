#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef __GENETIC__
#define __GENETIC__

#define POPSIZE 50//种群规模
#define MAXGENS 500//最大进化代数
#define NVARS 2//变量个数
#define PXOVER 0.25//交叉概率
#define PMUTATION 0.01//变异概率
#define TRUE 1
#define FALSE 0
#define PI 3.14159

/*
 * struct genotype
 * binary representation 
 *
 **/
struct genotype{
	char *gene;//代表染色体的基因位用0-1表示
	double fitness;//适应度函数
	double upper[NVARS];//最小值
	double lower[NVARS];//最大值
	int precision[NVARS];//准确率
	int p[NVARS];   //g 记录的是每段的基因段长度，即每个变量的编码长度
	double rfitness;//第一个目标函数
	double cfitness;//第二个目标函数
};

/*
 * Main Class
 * Standard Genetic Algorithm
 * To optimize real functions using binary representation
 * and standard operators
 * 
**/

class Genetic{
private:
	int GENESIZE;   //g 染色体编码长度
	int generation;  //g 迭代次数
	int cur_best;
	FILE *galog;  //g 过程记录文件
	char *filename;
	char *galog_file;
	int seed;  //随机数种子
	struct genotype population[POPSIZE+1];
	struct genotype newpopulation[POPSIZE+1];
public:
	Genetic(char *file, char* g_file, int i_seed): filename(file), galog_file(g_file), seed(i_seed) {};//构造函数，初始化必要参数：输入文件名称，过程记录文件名称，随机种子
	void initialize(void);
	inline double randval(double, double);
	void evaluate(void);
	void keep_the_best(void);
	void elitist(void);
	void select_proportional(void);
	void select_proportional_wheel(void);
	void select_tournament(void);
	void select_ranking(void);
	void crossover(void);
	void Xover(int, int);
	inline void swap(char*, char*);
	void mutate(void);
	void report(void);
	void solve(void);
	inline void bit_encode(double, double, double, int, char *);
	inline void bit_decode(char *, double, double, int, double &);
	inline int compute_p(double, double, int);
};
#endif
