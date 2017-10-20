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

void createNullMatrix(int** matr, int n) {
	for (int j = 0; j < n; j++) {
		matr[j] = new int[n];
		for (int k = 0; k < n; k++) {
			matr[j][k] = 0;
		}
	}
}

int** multMatrixParallelIn(int** matr1, int** matr2, int n, int threads_num) {
	int** matr = new int*[n];	createNullMatrix(matr, n);

	omp_set_num_threads(threads_num);
	double time1, time2;
	time1 = omp_get_wtime();
	int i = 0, j = 0;

#pragma omp parallel for shared(matr, matr1, matr2) private(i)
	for (i = 0; i < n; i++) {
#pragma omp parallel for shared(matr, matr1, matr2) private(j)
		for (j = 0; j < n; j++) for (int k = 0; k < n; k++) matr[i][j] += matr1[i][k] * matr2[k][j];	
	}
	time2 = omp_get_wtime();
	res_time = (time2 - time1) / (double)CLOCKS_PER_SEC;
	cout << "* Matrix mult.  (In)   Threads: " << threads_num << ", n = " << n << " time: " << time2 - time1 << " sec \n";
	return matr;
}

int** multMatrixParallelOut(int** matr1, int** matr2, int n, int threads_num) {
	int** matr = new int*[n]; createNullMatrix(matr, n);

	omp_set_num_threads(threads_num);
	double time1, time2; time1 = omp_get_wtime();
	int i = 0;
#pragma omp parallel for shared(matr, matr1, matr2) private(i)
	for (i = 0; i < n; i++) for (int j = 0; j < n; j++) for (int k = 0; k < n; k++) matr[i][j] += matr1[i][k] * matr2[k][j];
	
	time2 = omp_get_wtime(); res_time = (time2 - time1) / (double)CLOCKS_PER_SEC;
	cout << "* Matrix mult.  (Out)  Threads: " << threads_num << ", n = " << n << " time: " << time2 - time1 << " sec \n";
	return matr;
}



int** multmatrixParallelBlock(int** matr1, int** matr2, int n, int r, int threadsNum) {
	int** matr = new int*[n];	
	createNullMatrix(matr, n);
	int q = n / r;
	omp_set_num_threads(threadsNum);
	double time1, time2;
	time1 = omp_get_wtime();
	int i = 0, j = 0;
#pragma omp parallel for shared(matr, matr1, matr2) private(i)
	for (i = 0; i < q; i++) {
#pragma omp parallel 
		for (int j = 0; j < q; j++) 
			for (int k = 0; k < q; k++) 
				for (int i1 = i * r; i1 < (i + 1) * r; i1++) 
					for (int j1 = j * r; j1 < (j + 1) * r; j1++) 
						for (int k1 = k * r; k1 < (k + 1) * r; k1++) 
							matr[i1][j1] += matr1[i1][k1] * matr2[k1][j1];		
	}
	time2 = omp_get_wtime();
	cout << "* Matrix mult. (block) Threads: " << threadsNum  << ", n = " << n << " time: " << time2 - time1 << " sec, r = " << r << ", block number (Q) = " << q << endl;
	return matr;
}

void generateMatrix(int** matr1, int** matr2, int n) {
	int first_value = -100;
	int last_value = 200;
	for (int i = 0; i < n; i++) {
		matr1[i] = new int[n];
		matr2[i] = new int[n];
		for (int j = 0; j < n; j++) {
			matr1[i][j] = first_value + rand() % last_value;
			matr2[i][j] = first_value + rand() % last_value;
		}
	}
}

void main(int argc, char** argv) {
	int block_num = 100;
	int thread_num = 10;
	
	int block_nums[5] = { 10, 50, 100, 250, 500 };

	for (int n = 100; n < 5000; n += 500) {
		cout << "Matrix size: " << n << endl;
		int** matr1 = new int*[n];
		int** matr2 = new int*[n];
		int** res_matr1, **res_matr2;
		
		generateMatrix(matr1, matr2, n);

		for (thread_num = 1; thread_num < n / 10; thread_num += 10) {
			
			multMatrixParallelIn(matr1, matr2, n, thread_num);
			multMatrixParallelOut(matr1, matr2, n, thread_num);
			for (int r = 1; r < n/2; r += 50) {
				multmatrixParallelBlock(matr1, matr2, n, r, thread_num);
			}
		}
		
	}
	system("pause");

}











































/*int** matrix_multipl(int** matr1, int** matr2, int n) {
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
}*/


/*int** matrix_multipl_par_two(int** matr1, int** matr2, int n, int threads_num) {
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
}*/








































/*// L1.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include <iostream>
#include <cstdlib>
#include <omp.h>

using namespace std;

void randomiseMatrix(int **matrix, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] = rand() % 11;
		}
	}

	return;
}

void printMatrix(int **matrix, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << matrix[i][j] << "\t";
		}
		cout << "\n";
	}
	cout << "\n\n";
	return;
}


void pointCalc(int **matrix1, int **matrix2, int **result, int n, int threadsNum) {
	//Устанавливаем число потоков

	omp_set_num_threads(threadsNum);
	int i, j, k;

#pragma omp parallel for shared(matrix1, matrix2, result) private(i, j, k)
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			result[i][j] = 0;
			for (k = 0; k < n; k++) {
				result[i][j] += (matrix1[i][k] * matrix2[k][j]);
			}
		}
	}
	return;
}

void Rec_Mult(int *C, const int *A, const int *B, int n, int rowsize);



void multRec(int** res, int** a, int** b, int n, int r) {
	int q = n / r;
	for (int i = 0; i < q-1; i++) {
		for (int j = 0; j < q-1; j++) {
			//res[i][j] = 0;
			for (int k = 0; k < q-1; k++) {
				res[i+1][j+1] += (a[i+1][k+1] * b[k+1][j+1]);
			}
		}
	}
}

int main(int argc, char** argv) {
	//srand(time(NULL));
	cout << "\nEnter size of matrices: "; int n = 1000; cin >> n;

	int threadsNum = 2;
	cout << "\nEnter number of threads: "; cin >> threadsNum;
	

	int **matrix1;
	int **matrix2;

	matrix1 = (int**)malloc(sizeof(int*)*n);
	for (int i = 0; i < n; i++) {
		matrix1[i] = (int*)malloc(sizeof(int)*n);
	}
	matrix2 = (int**)malloc(sizeof(int*)*n);
	for (int i = 0; i < n; i++) {
		matrix2[i] = (int*)malloc(sizeof(int)*n);
	}

	//Генерируем случайные матрицы для умножения
	randomiseMatrix(matrix1, n);
	randomiseMatrix(matrix2, n);

	cout << "\nMatrix1:\n";
	printMatrix(matrix1, n);
	cout << "\nMatrix2:\n";
	printMatrix(matrix2, n);

	int **result = (int**)malloc(sizeof(int*)*n);;
	for (int i = 0; i < n; i++) {
		result[i] = (int*)malloc(sizeof(int)*n);
	}

	pointCalc(matrix1, matrix2, result, n, threadsNum);

	int *c = (int*)malloc(sizeof(int)*n);
	int *a = (int*)malloc(sizeof(int)*n);
	int *b = (int*)malloc(sizeof(int)*n);

	for (int z = 0, i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a[z] = matrix1[i][j]; b[z] = matrix2[i][j]; c[z] = 0; z++;
		}
	}

	

	cout << "\nResult:\n";
	printMatrix(result, n);


	int **result2 = (int**)malloc(sizeof(int*)*n);;
	for (int i = 0; i < n; i++) {
		result2[i] = (int*)malloc(sizeof(int)*n);
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			result2[i][j] = 0;
		}
	}
	
	multRec(result2, matrix1, matrix2, n, 2);


	cout << "\nResult(block):\n";
	printMatrix(result2, n);

	/*
	Rec_Mult(c, a, b, n, n);

	int **result2 = (int**)malloc(sizeof(int*)*n);;
	for (int i = 0; i < n; i++) {
		result2[i] = (int*)malloc(sizeof(int)*n);
	}
	for (int z = 0, i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			result2[i][j] = c[z];
			z++;
		}
	}
	cout << "\nResult(block):\n";
	printMatrix(result2, n);
	*//*
	return 0;
}


void Rec_Mult(int *C, const int *A, const int *B, int n, int rowsize)
{
	if (n == 2)
	{
		const int d11 = 0;
		const int d12 = 1;
		const int d21 = rowsize;
		const int d22 = rowsize + 1;

		C[d11] += A[d11] * B[d11] + A[d12] * B[d21];
		C[d12] += A[d11] * B[d12] + A[d12] * B[d22];
		C[d21] += A[d21] * B[d11] + A[d22] * B[d21];
		C[d22] += A[d21] * B[d12] + A[d22] * B[d22];
	}
	else
	{
		const int d11 = 0;
		const int d12 = n / 2;
		const int d21 = (n / 2) * rowsize;
		const int d22 = (n / 2) * (rowsize + 1);

		// C11 += A11 * B11
		Rec_Mult(C + d11, A + d11, B + d11, n / 2, rowsize);
		// C11 += A12 * B21
		Rec_Mult(C + d11, A + d12, B + d21, n / 2, rowsize);

		// C12 += A11 * B12
		Rec_Mult(C + d12, A + d11, B + d12, n / 2, rowsize);
		// C12 += A12 * B22
		Rec_Mult(C + d12, A + d12, B + d22, n / 2, rowsize);

		// C21 += A21 * B11
		Rec_Mult(C + d21, A + d21, B + d11, n / 2, rowsize);
		// C21 += A22 * B21
		Rec_Mult(C + d21, A + d22, B + d21, n / 2, rowsize);

		// C22 += A21 * B12
		Rec_Mult(C + d22, A + d21, B + d12, n / 2, rowsize);
		// C22 += A22 * B22
		Rec_Mult(C + d22, A + d22, B + d22, n / 2, rowsize);
	}
}

void Rec_Mult_R(double* C, int offC, double* A, int offA, double* B, int offB, int n, int rowsize) {
	if (n == 1) {
		C[offC] += A[offA] * B[offB];
	}
	else {
		int d11 = 0;
		int d12 = n / 2;
		int d21 = (n / 2) * rowsize;
		int d22 = (n / 2) * (rowsize + 1);

		int C11 = offC + d11;
		int A11 = offA + d11;
		int B11 = offB + d11;

		int C12 = offC + d12;
		int A12 = offA + d12;
		int B12 = offB + d12;

		int C21 = offC + d21;
		int A21 = offA + d21;
		int B21 = offB + d21;

		int C22 = offC + d22;
		int A22 = offA + d22;
		int B22 = offB + d22;

		// C11 += A11 * B11
		Rec_Mult_R(C, C11, A, A11, B, B11, n / 2, rowsize);
		// C11 += A12 * B21
		Rec_Mult_R(C, C11, A, A12, B, B21, n / 2, rowsize);

		// C12 += A11 * B12
		Rec_Mult_R(C, C12, A, A11, B, B12, n / 2, rowsize);
		// C12 += A12 * B22
		Rec_Mult_R(C, C12, A, A12, B, B22, n / 2, rowsize);

		// C21 += A21 * B11
		Rec_Mult_R(C, C21, A, A21, B, B11, n / 2, rowsize);
		// C21 += A22 * B21
		Rec_Mult_R(C, C21, A, A22, B, B21, n / 2, rowsize);

		// C22 += A21 * B12
		Rec_Mult_R(C, C22, A, A21, B, B12, n / 2, rowsize);
		// C22 += A22 * B22
		Rec_Mult_R(C, C22, A, A22, B, B22, n / 2, rowsize);
	}
}*/