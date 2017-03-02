/*
���ѡ�񣬵��㽻�棬������죬��������Ϊ����С�����ֵ����ֵ��0.01
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#define upper_bound 10 //�����������Ͻ�
#define lower_bound 1 //�����������½�
#define popu_size 400 //��Ⱥ��С
#define chromosome_size 17 //Ⱦɫ�峤�ȣ����ݾ���ΪС�������λ�����
#define gene_max_times (long) 2 << 16//����������
#define survival_rate 0.6 //����ʣ�����ѡ�����Ⱥ�ĸ�����
#define mutate_rate 0.1 //�������
int gene_times = 0;//��������
int survival_amount = popu_size * survival_rate;//��Ⱥ����ѡ���ʣ�������
int popu[popu_size][chromosome_size] = { 0 };//��Ⱥ��ÿһ�б�ʾһ������,ȫ����ʼ��Ϊ0
int new_popu[popu_size][chromosome_size] = { 0 };//ѡ������в�������һ����Ⱥ
int current_best[chromosome_size] = { 0 };//���ڴ�ŵ�ǰ��Ⱥ���Ž��
double fitness[popu_size] = { 0 };//������Ӧ��,�ڼ�����Сֵʱ��fitness����Ԫ��Ϊ����ֵ���෴��

/*�˺������ڳ�ʼ����Ⱥ�и������Ⱦɫ������*/
void init_popu() {
	
	for (int i = 0; i < popu_size; i++) {
		for (int j = 0; j < chromosome_size; j++) {
			if (rand() > RAND_MAX/2) {
				popu[i][j] = 1;
			}
		}
	}
}

/*�˺�����Ⱦɫ���Ӧ����ת��Ϊ�����������ж�Ӧ��ֵ*/
double chromosome2num(int *chro) {
	int to_return = 0;
	for (int i = 0; i < chromosome_size ; i++) {
		to_return = to_return << 1;
		if (chro[i]) {
			to_return = to_return + 1;
		}
	}//��ʱ�ѽ�Ⱦɫ���Ӧ����ת��Ϊһ��int������
	return (lower_bound + (upper_bound - lower_bound) * to_return * 1.0 / (2 << chromosome_size - 1));//�˴�����1.0��Ϊ��������Ϊ������
}

/*��ǰҪ�����ֵ�ĺ���*/
double func(double value) {
	return value * cos(value);
}

/*�˺���������Ⱥ�и��������Ӧ�ȣ��ڼ�����Сֵʱ��fitness����Ԫ��Ϊ����ֵ���෴��*/
void cal_fitness(int direction) {
	for (int i = 0; i < popu_size; i++) {
		double current = chromosome2num(popu[i]);
		fitness[i] = func(current);
		if (direction) {
			fitness[i] = -fitness[i];
		}
	}
}

/*�˺������ڱ�����Ƹ���*/
void mem_cpy(int* dest, int* src) {
	for (int i = 0; i < chromosome_size; i++) {
		dest[i] = src[i];
	}
}

/*�˺�������������ʱ������Ⱦɫ��λ�û����ڽ���ʱ�����彻������Ⱦɫ�壬����Ϊ�������������ݣ�cursorΪȾɫ�彻������ʼλ��*/
void mem_exchange(int*dest, int* src, int cursor, int ends) {
	int tmp_for_exchange = 0;
	for (int i = cursor; i < ends; i++) {
		tmp_for_exchange = dest[i];
		dest[i] = src[i];
		src[i] = tmp_for_exchange;
	}
}

/*�˺���������Ⱥ����Ⱦɫ�彻��λ��*/
void chromosome_exchange(int i, int j) {
	double tmp_for_fitness_change = 0.0;
	mem_exchange(popu[i], popu[j], 0, chromosome_size);
	tmp_for_fitness_change = fitness[i];
	fitness[i] = fitness[j];
	fitness[j] = tmp_for_fitness_change;
}

/*�˺����Ը�������Ӧ�Ƚ�������,����������Ϊfalse���򰴴�С�������򣬷��򰴴Ӵ�С����*/
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

/*�˺���������Ч�ı�ţ�����������һ��0��up_limit֮��������*/
int get_number(int up_limit) {
	return rand() * 1.0 * (up_limit - 1) / RAND_MAX;
}

/*���ѡ������Բ�����һ����Ⱥ��ѡ������Ϊ�������Դ���ʣ�ʣ���λ���ɱ�ѡ��ĸ��彻�����,����ֵΪѡ�����Ⱥ����*/
void select() {
	for (int i = 0; i < survival_amount; i++) {
		int to_select = get_number(popu_size);
		while (to_select < popu_size * 0.1) {
			to_select = get_number(popu_size);
		}
		mem_cpy(new_popu[i], popu[to_select]);
	}
}

/*��ԭ��Ⱥ�����ѡ�������н���*/
void cross_over() {
	int left_seat = popu_size - survival_amount;//��̭�ĸ������������Ҫ�������ɵĸ��������
	int count_to = (left_seat % 2 == 0) ? popu_size : popu_size - 1;//Ϊ���ø���ɶԽ��棬�Խ����������΢��
	int is_changed[popu_size] = { 0 };//���ڼ�¼ԭ��Ⱥ�ж�Ӧ�����Ƿ��ѽ����
	for (int i = survival_amount; i < count_to; i = i + 2) {
		int father = get_number(popu_size);
		int mother = get_number(popu_size);
		while (father == mother || is_changed[father] || is_changed[mother]) {//��֤ÿ�ν����˫����ͬ���Ҷ��ǵ�һ�ν���
			if (is_changed[father]) {
				father = get_number(popu_size);
			}
			if (is_changed[mother] || father == mother) {
				mother = get_number(popu_size);
			}
		}
		int cross_after = get_number(chromosome_size);
		while (cross_after < 1 || cross_after > chromosome_size - 2) {//��֤����Ĳ��ֲ���̫���̫�٣���Ҫ��ֹ��ȫ��������ȫ������
			cross_after = get_number(chromosome_size);
		}
		//���н���,�����Ƶ�����Ⱥ������
		mem_exchange(popu[father], popu[mother], cross_after, chromosome_size);
		mem_cpy(new_popu[i], popu[father]);
		mem_cpy(new_popu[i + 1], popu[mother]);
		is_changed[father] = 1;
		is_changed[mother] = 1;
	}
}

/*���죬����Ⱥ��ÿ�����嶼��һ�����ʱ���*/
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

/*������Ⱥ�еĸ��帴�Ƶ�ԭ��Ⱥ��*/
void cpy_from_new_to_old() {
	for (int i = 0; i < popu_size - 1; i++) {
		mem_cpy(popu[i], new_popu[i]);
	}
}

/*�ж��Ƿ�����*/
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

/*��ѭ��,����Ϊ0��Ϊ�����ֵ������Ϊ1��Ϊ����Сֵ*/
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
	printf("����㣺%f, ����ֵ��%f, ��������:%d\n", max, max_value, gene_times);
	double min = gene_algorithm(1);
	double min_value = func(min);
	printf("��С�㣺%f, ����ֵ��%f, ��������:%d", min, min_value, gene_times);
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