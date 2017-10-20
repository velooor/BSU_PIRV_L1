#include "stdafx.h"

#include "stdlib.h"
#include <random>
#include <iostream>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <iomanip>

using namespace std;

double res_time;

int** matrix_multipl(int** matr1, int** matr2, int n) {
	int** matr = new int*[n];
	for (int j = 0; j < n; j++) {
		matr[j] = new int[n];
		for (int k = 0; k < n; k++) {
			matr[j][k] = 0;
		}
	}
	clock_t time1, time2;
	time1 = clock();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				matr[i][j] += matr1[i][k] * matr2[k][j];
			}
		}
	}
	time2 = clock();
	res_time = (time2 - time1) / (double)CLOCKS_PER_SEC;
	cout << "Time standart: " << (time2 - time1) / (double)CLOCKS_PER_SEC << endl;
	return matr;
}


int** matrix_multipl_block(int** matr1, int** matr2, int n, int block_size) {
	int** matr = new int*[n];
	int block_num = n / block_size;
	for (int j = 0; j < n; j++) {
		matr[j] = new int[n];
		for (int k = 0; k < n; k++) {
			matr[j][k] = 0;
		}
	}
	clock_t time1, time2;
	time1 = clock();
	for (int i = 0; i < block_num; i++) {
		for (int j = 0; j < block_num; j++) {
			for (int k = 0; k < block_num; k++) {
				for (int i1 = i * block_size; i1 < (i + 1) * block_size; i1++) {
					for (int j1 = j * block_size; j1 < (j + 1) * block_size; j1++) {
						for (int k1 = k * block_size; k1 < (k + 1) * block_size; k1++) {
							matr[i1][j1] += matr1[i1][k1] * matr2[k1][j1];
						}
					}
				}
			}
		}
	}
	time2 = clock();
	res_time = (time2 - time1) / (double)CLOCKS_PER_SEC;
	cout << "Time standart block " << block_num << ": " << (time2 - time1) / (double)CLOCKS_PER_SEC << endl;
	return matr;
}

int** matrix_multipl_par(int** matr1, int** matr2, int n, int threads_num) {
	int** matr = new int*[n];
	for (int j = 0; j < n; j++) {
		matr[j] = new int[n];
		for (int k = 0; k < n; k++) {
			matr[j][k] = 0;
		}
	}
	omp_set_num_threads(threads_num);
	double time1, time2;
	time1 = omp_get_wtime();
	int i = 0;
#pragma omp parallel for shared(matr, matr1, matr2) private(i)
	for (i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				matr[i][j] += matr1[i][k] * matr2[k][j];
			}
		}
	}
	time2 = omp_get_wtime();
	res_time = (time2 - time1) / (double)CLOCKS_PER_SEC;
	cout << "Time paral(par) " << threads_num << ":" << time2 - time1 << '\n';
	return matr;
}

int** matrix_multipl_par_two(int** matr1, int** matr2, int n, int threads_num) {
	int** matr = new int*[n];
	for (int j = 0; j < n; j++) {
		matr[j] = new int[n];
		for (int k = 0; k < n; k++) {
			matr[j][k] = 0;
		}
	}
	omp_set_num_threads(threads_num);
	double time1, time2;
	time1 = omp_get_wtime();
	int i = 0;
	int j = 0;
#pragma omp parallel for shared(matr, matr1, matr2) private(i)
	for (i = 0; i < n; i++) {
#pragma omp parallel for shared(matr, matr1, matr2) private(j)
		for (j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				matr[i][j] += matr1[i][k] * matr2[k][j];
			}
		}
	}
	time2 = omp_get_wtime();
	res_time = (time2 - time1) / (double)CLOCKS_PER_SEC;
	cout << "Time paral(par_two) " << threads_num << ":" << time2 - time1 << '\n';
	return matr;
}

int** matrix_multipl_block_par(int** matr1, int** matr2, int n, int block_size, int threads_num) {
	int** matr = new int*[n];
	int block_num = n / block_size;
	for (int j = 0; j < n; j++) {
		matr[j] = new int[n];
		for (int k = 0; k < n; k++) {
			matr[j][k] = 0;
		}
	}
	omp_set_num_threads(threads_num);
	double time1, time2;
	time1 = omp_get_wtime();
	int i = 0;
	int j = 0;
#pragma omp parallel for shared(matr, matr1, matr2) private(i)
	for (i = 0; i < block_num; i++) {
#pragma omp parallel 
		for (int j = 0; j < block_num; j++) {
			for (int k = 0; k < block_num; k++) {
				for (int i1 = i * block_size; i1 < (i + 1) * block_size; i1++) {
					for (int j1 = j * block_size; j1 < (j + 1) * block_size; j1++) {
						for (int k1 = k * block_size; k1 < (k + 1) * block_size; k1++) {
							matr[i1][j1] += matr1[i1][k1] * matr2[k1][j1];
						}
					}
				}
			}
		}
	}
	time2 = omp_get_wtime();
	cout << "Time paral(block) " << threads_num << " block_num - " << block_num << ":" << time2 - time1 << '\n';
	return matr;
}

void main(int argc, char** argv) {
	int n = 500;
	if (argv[1] != NULL)
		n = *argv[1];
	int block_num = 100;
	int thread_num = 10;
	int firs_value = -100;
	int last_value = 200;
	int block_nums[5] = { 10, 50, 100, 250, 500 };
	ofstream myfile;
	myfile.open("res.csv");
	myfile << "n,block_num,thread_num,time\n";
	for (n = 500; n < 3000; n += 500) {
		cout << n << endl;
		int** matr1 = new int*[n];
		int** matr2 = new int*[n];
		int** res_matr1, **res_matr2;
		for (int i = 0; i < n; i++) {
			matr1[i] = new int[n];
			matr2[i] = new int[n];
			for (int j = 0; j < n; j++) {
				matr1[i][j] = firs_value + rand() % last_value;
				matr2[i][j] = firs_value + rand() % last_value;
			}
		}
		res_matr1 = matrix_multipl(matr1, matr2, n);
		myfile << n << "," << 1 << "," << 1 << "," << fixed << setprecision(5) << res_time << endl;
		for (thread_num = 2; thread_num < n / 10; thread_num += 20) {
			res_matr1 = matrix_multipl_par(matr1, matr2, n, thread_num);
			myfile << n << "," << 1 << "," << thread_num << "," << fixed << setprecision(5) << res_time << endl;
			res_matr1 = matrix_multipl_par_two(matr1, matr2, n, thread_num);
			myfile << n << "," << 1 << "," << thread_num << "," << fixed << setprecision(5) << res_time << endl;
			for (int k = 0; k < 5; k++) {
				block_num = block_nums[k];
				res_matr1 = matrix_multipl_block_par(matr1, matr2, n, n / block_num, thread_num);
				myfile << n << "," << block_num << "," << thread_num << "," << fixed << setprecision(5) << res_time << endl;
			}
			if (thread_num == 2)
				thread_num = 0;
		}
		for (int k = 0; k < 5; k++) {
			block_num = block_nums[k];
			res_matr2 = matrix_multipl_block(matr1, matr2, n, n / block_num);
			myfile << n << "," << block_num << "," << 1 << "," << fixed << setprecision(5) << res_time << endl;
		}
	}
	myfile.close();
	system("pause");
}