#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ran2.c"

#define eqs_MAX 1000
#define mcs_MAX 1000000
#define L_MAX 32

long seed = -1;

void init_lattice(int L, int N, int** lattice, int** idx); // 초기화 (랜덤)
void random_index(int L, int N, int* rand_idx); // 랜덤 인덱스 리스트 생성
double compute_energy(int L, int N, int** lattice, int i, int j); // site 에너지 계산
void monte_carlo(int L, int N, int** lattice, int** idx, int* rand_idx, FILE* file); // 몬테 카를로 방법으로 observables 계산

int main() {
	FILE* file;
	int i, L, N, ** lattice, ** idx, * rand_idx;
	double T;

	// observables 파일 열기
	file = fopen("U.txt", "w");
	if (file == NULL) {
		printf("file ERROR\n");
		exit(1);
	}
	fprintf(file, "L\t");

	for (T = 2.3; T > 2.2; T -= 0.002) {
		fprintf(file, "%f\t", T);
	}
	fprintf(file, "\n");

	// monte_carlo
	clock_t t0 = clock(); // 반복문 시작 시간

	for (L = 4; L <= L_MAX; L *= 2) {

		clock_t t1 = clock(); // 시작 시간

		N = L * L;

		// 동적할당
		lattice = (int**)malloc(sizeof(int*) * L);
		for (i = 0; i < L; i++) {
			lattice[i] = (int*)malloc(sizeof(int) * L);
		}
		idx = (int**)malloc(sizeof(int*) * N);
		for (i = 0; i < N; i++) {
			idx[i] = (int*)malloc(sizeof(int) * 2);
		}
		rand_idx = (int*)malloc(sizeof(int) * N);

		monte_carlo(L, N, lattice, idx, rand_idx, file);

		// 동적할당 해제
		for (i = 0; i < L; i++) {
			free(lattice[i]);
		}
		free(lattice);
		free(idx);
		free(rand_idx);

		clock_t t1_end = clock();
		printf("\t%.3f s\n", (double)(t1_end - t1) / CLOCKS_PER_SEC); // 소요 시간
	}
	clock_t t0_end = clock();
	printf("total elapsed time : %.3f s\n", (double)(t0_end - t0) / CLOCKS_PER_SEC); // 반복문 소요 시간

	fclose(file);

	return 0;
}

void init_lattice(int L, int N, int** lattice, int** idx) { // 초기화 (랜덤)
	int i, j, k;

	k = 0;
	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			if (ran2(&seed) > 0.5) {
				lattice[i][j] = 1;
			}
			else lattice[i][j] = -1;

			idx[k][0] = i;
			idx[k][1] = j;

			k++;
		}
	}
}

void random_index(int L, int N, int* rand_idx) { // 랜덤 인덱스 리스트 생성
	int cnt, num, * tmp;

	tmp = (int*)calloc(sizeof(int), N);

	cnt = 0;
	while (cnt < N) {

		num = (int)(ran2(&seed) * N);
		if (tmp[num] == 0) {
			rand_idx[cnt] = num;
			tmp[num] = 1;
			cnt++;
		}
	}

	free(tmp);
}

double compute_energy(int L, int N, int** lattice, int i, int j) { // site 에너지 계산
	int LEFT, RIGHT, ABOVE, BELOW;

	// neighbors & periodic boundary condition
	LEFT = lattice[(i - 1 + L) % L][j];
	RIGHT = lattice[(i + 1) % L][j];
	ABOVE = lattice[i][(j - 1 + L) % L];
	BELOW = lattice[i][(j + 1) % L];

	return -lattice[i][j] * (LEFT + RIGHT + ABOVE + BELOW);
}

void monte_carlo(int L, int N, int** lattice, int** idx, int* rand_idx, FILE* file) { // 몬테 카를로 방법으로 |M| 계산
	int i, j, k, eqs, mcs;
	double T, M, sum_M2, sum_M4, E1, E0, dE, U;

	fprintf(file, "%d\t", L);
	printf("%d\t", L);

	init_lattice(L, N, lattice, idx);

	// init_M
	M = 0;
	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			M += lattice[i][j];
		}
	}

	// monte_carlo
	for (T = 2.3; T > 2.2; T -= 0.002) {
		// eqs
		for (eqs = 0; eqs < eqs_MAX; eqs++) {

			random_index(L, N, rand_idx);
			for (k = 0; k < N; k++) {
				i = idx[rand_idx[k]][0];
				j = idx[rand_idx[k]][1];

				E0 = compute_energy(L, N, lattice, i, j);
				lattice[i][j] *= -1;
				E1 = compute_energy(L, N, lattice, i, j);

				dE = E1 - E0;

				if (dE > 0 && ran2(&seed) > exp(-dE / T)) lattice[i][j] *= -1; // skip
				else M += (double)lattice[i][j] * 2; // flip
			}
		}

		// mcs
		sum_M2 = 0;
		sum_M4 = 0;
		for (mcs = 0; mcs < mcs_MAX; mcs++) {

			random_index(L, N, rand_idx);
			for (k = 0; k < N; k++) {
				i = idx[rand_idx[k]][0];
				j = idx[rand_idx[k]][1];

				E0 = compute_energy(L, N, lattice, i, j);
				lattice[i][j] *= -1;
				E1 = compute_energy(L, N, lattice, i, j);

				dE = E1 - E0;

				if (dE > 0 && ran2(&seed) > exp(-dE / T)) lattice[i][j] *= -1; // skip
				else M += (double)lattice[i][j] * 2; // flip
			}
			sum_M2 += (M * M);
			sum_M4 += (M * M * M * M);
		}
		U = 1 - ((sum_M4 / (double)mcs_MAX) / (3 * (sum_M2 / (double)mcs_MAX) * (sum_M2 / (double)mcs_MAX)));
		fprintf(file, "%f\t", U);
		printf("*");
	}
	fprintf(file, "\n");
}
