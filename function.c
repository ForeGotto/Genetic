/*
随机选择，单点交叉，随机变异，收敛条件为极差小于最大值绝对值的0.01
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#define upper_bound 10 //定义域区间上界
#define lower_bound 1 //定义域区间下界
#define popu_size 400 //种群大小
#define chromosome_size 17 //染色体长度，根据精度为小数点后四位计算得
#define gene_max_times (long) 2 << 16//最大迭代次数
#define survival_rate 0.6 //存活率，决定选择后种群的个体数
#define mutate_rate 0.1 //变异概率
int gene_times = 0;//迭代次数
int survival_amount = popu_size * survival_rate;//种群进行选择后剩余的数量
int popu[popu_size][chromosome_size] = { 0 };//种群，每一行表示一个个体,全部初始化为0
int new_popu[popu_size][chromosome_size] = { 0 };//选择过程中产生的下一代种群
int current_best[chromosome_size] = { 0 };//用于存放当前种群最优结果
double fitness[popu_size] = { 0 };//个体适应度,在计算最小值时，fitness数组元素为函数值的相反数

/*此函数用于初始化种群中各个体的染色体数组*/
void init_popu() {
	
	for (int i = 0; i < popu_size; i++) {
		for (int j = 0; j < chromosome_size; j++) {
			if (rand() > RAND_MAX/2) {
				popu[i][j] = 1;
			}
		}
	}
}

/*此函数把染色体对应数组转换为函数定义域中对应的值*/
double chromosome2num(int *chro) {
	int to_return = 0;
	for (int i = 0; i < chromosome_size ; i++) {
		to_return = to_return << 1;
		if (chro[i]) {
			to_return = to_return + 1;
		}
	}//此时已将染色体对应数组转化为一个int型整数
	return (lower_bound + (upper_bound - lower_bound) * to_return * 1.0 / (2 << chromosome_size - 1));//此处乘以1.0是为了运算结果为浮点数
}

/*当前要求最大值的函数*/
double func(double value) {
	return value * cos(value);
}

/*此函数计算种群中各个体的适应度，在计算最小值时，fitness数组元素为函数值的相反数*/
void cal_fitness(int direction) {
	for (int i = 0; i < popu_size; i++) {
		double current = chromosome2num(popu[i]);
		fitness[i] = func(current);
		if (direction) {
			fitness[i] = -fitness[i];
		}
	}
}

/*此函数用于变异后复制个体*/
void mem_cpy(int* dest, int* src) {
	for (int i = 0; i < chromosome_size; i++) {
		dest[i] = src[i];
	}
}

/*此函数用于在排序时交换两染色体位置或者在交叉时两个体交换部分染色体，功能为交换两数组内容，cursor为染色体交换的起始位置*/
void mem_exchange(int*dest, int* src, int cursor, int ends) {
	int tmp_for_exchange = 0;
	for (int i = cursor; i < ends; i++) {
		tmp_for_exchange = dest[i];
		dest[i] = src[i];
		src[i] = tmp_for_exchange;
	}
}

/*此函数用于种群中两染色体交换位置*/
void chromosome_exchange(int i, int j) {
	double tmp_for_fitness_change = 0.0;
	mem_exchange(popu[i], popu[j], 0, chromosome_size);
	tmp_for_fitness_change = fitness[i];
	fitness[i] = fitness[j];
	fitness[j] = tmp_for_fitness_change;
}

/*此函数对各个体适应度进行排序,如果传入参数为false，则按从小到大排序，否则按从大到小排序*/
void rank_popu() {
		for (int i = 0; i < popu_size; i++) {
		for (int j = i + 1; j < popu_size; j++) {
			if (fitness[i] > fitness[j]) {
				chromosome_exchange(i, j);
			}
		}
	}
		mem_cpy(current_best, popu[popu_size - 1]);
}

/*此函数返回有效的编号，功能是生成一个0到up_limit之间的随机数*/
int get_number(int up_limit) {
	return rand() * 1.0 * (up_limit - 1) / RAND_MAX;
}

/*随机选择个体以产生下一代种群，选择数量为总数乘以存活率，剩余的位置由被选择的个体交叉产生,返回值为选择后种群数量*/
void select() {
	for (int i = 0; i < survival_amount; i++) {
		int to_select = get_number(popu_size);
		while (to_select < popu_size * 0.1) {
			to_select = get_number(popu_size);
		}
		mem_cpy(new_popu[i], popu[to_select]);
	}
}

/*从原种群中随机选择个体进行交叉*/
void cross_over() {
	int left_seat = popu_size - survival_amount;//淘汰的个体的数量，即要交叉生成的个体的数量
	int count_to = (left_seat % 2 == 0) ? popu_size : popu_size - 1;//为了让个体成对交叉，对交叉对数进行微调
	int is_changed[popu_size] = { 0 };//用于记录原种群中对应个体是否已交叉过
	for (int i = survival_amount; i < count_to; i = i + 2) {
		int father = get_number(popu_size);
		int mother = get_number(popu_size);
		while (father == mother || is_changed[father] || is_changed[mother]) {//保证每次交叉的双方不同，且都是第一次交叉
			if (is_changed[father]) {
				father = get_number(popu_size);
			}
			if (is_changed[mother] || father == mother) {
				mother = get_number(popu_size);
			}
		}
		int cross_after = get_number(chromosome_size);
		while (cross_after < 1 || cross_after > chromosome_size - 2) {//保证交叉的部分不会太多或太少，主要防止完全交换或完全不交叉
			cross_after = get_number(chromosome_size);
		}
		//进行交叉,并复制到新种群数组中
		mem_exchange(popu[father], popu[mother], cross_after, chromosome_size);
		mem_cpy(new_popu[i], popu[father]);
		mem_cpy(new_popu[i + 1], popu[mother]);
		is_changed[father] = 1;
		is_changed[mother] = 1;
	}
}

/*变异，新种群中每个个体都有一定几率变异*/
void mutate() {
	for (int i = 0; i < popu_size - 2; i++) {
		double current_rate = rand() * 1.0 / RAND_MAX;
		if (current_rate <= mutate_rate) {
			int mutate_position = get_number(chromosome_size);
			
			while (mutate_position < 2 || mutate_position > chromosome_size - 2) {
				mutate_position = get_number(chromosome_size);
			}
			new_popu[i][mutate_position] = 1 - new_popu[i][mutate_position];
		}
	}
}

/*把新种群中的个体复制到原种群中*/
void cpy_from_new_to_old() {
	for (int i = 0; i < popu_size - 1; i++) {
		mem_cpy(popu[i], new_popu[i]);
	}
}

/*判断是否收敛*/
int is_converged() {
	if ((fitness[popu_size - 1] - fitness[0]) < 0.01) {
		return 1;
	}
	else {
		return 0;
	}
}

void print_all() {
	printf("---------------\n");
	for (int i = 0; i < popu_size; i++) {
		for (int j = 0; j < chromosome_size; j++) {
			printf("%d", popu[i][j]);
		}
		double current = chromosome2num(popu[i]);
		printf(":%f,%f\n", current, fitness[i]);
	}
}

/*主循环,参数为0则为求最大值，参数为1则为求最小值*/
double gene_algorithm(int mode) {
	gene_times = 0;
	init_popu();
	cal_fitness(mode);
	rank_popu();
	double first_max = mode ? -fitness[popu_size - 1] : fitness[popu_size - 1];
	printf("%f\n", first_max);
	while (!is_converged() && gene_times < gene_max_times) {
	select();
	cross_over();
	mutate();
	cpy_from_new_to_old();
	cal_fitness(mode);
	rank_popu();
	gene_times++;
	//printf("%f\n", fitness[chromosome_size]);
	}
	double to_return = chromosome2num(popu[popu_size - 1]);
	return to_return;
}



void print_new() {
	for (int i = 0; i < popu_size; i++) {
		for (int j = 0; j < chromosome_size; j++) {
			printf("%d", new_popu[i][j]);
		}
		double current = chromosome2num(new_popu[i]);
		printf(":%f\n", current);
	}
}

int main() {
	srand(time(NULL));
	double max =  gene_algorithm(0);
	double max_value = func(max);
	printf("极大点：%f, 函数值：%f, 迭代次数:%d\n", max, max_value, gene_times);
	double min = gene_algorithm(1);
	double min_value = func(min);
	printf("极小点：%f, 函数值：%f, 迭代次数:%d", min, min_value, gene_times);
	/**
	//print_all();
	init_popu();
	cal_fitness(0);
	//print_all();
	rank_popu();

	//print_all();

	select();
	cross_over();
	mutate();
	cpy_from_new_to_old();
	cal_fitness(0);
	rank_popu();
	print_all();
	print_new();
	**/
	system("pause");
}