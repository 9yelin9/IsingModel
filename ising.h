// ising.h : headers for 2D ising model

#ifndef ISING_H
#define ISING_H

#define USE_MATH_DEFINES
#define OMP_THREAD (1)

#define T_MAX (51) // max temperature
#define EQS (1000) // equilibrium steps
#define MCS (10000) // monte carlo steps

#define ENERGY (s.e)
#define HEAT_CAPACITY(T) ((s.e2 - pow(s.e, 2)) / pow(T, 2))
#define ABS_MAGNETIZATION (s.ma)
#define MAG_SUSCEPTIBILITY(T) ((s.m2 - pow(s.ma, 2)) / T)
#define CUMULANT (1 - (s.m4 / (3 * s.m2)))

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct MonteCarloSolution {
	double m;  // magnetization
	double ma; // absolute magnetization
	double m2; // magnetization^2
	double m4; // magnetization^4
	double e;  // energy
	double e2; // energy^2
} Solution;

FILE* OpenFile(char *fs, char *ftype); // open file
void GenIdx(int N, int *idx); // generate random index array
void GenName(int L, char *data_type, char *fs); // generate random index array

void InitLatRand(int L, int (*lat)[L]); // initialize lattice randomly
void InitLatBin(int dec, int L, int (*lat)[L]); // initialize lattice by converting decimal to binary
void InitMonteCarlo(int L, int N, double T, int (*lat)[L]); // initialize Monte Carlo to reach a equilibrium state

double CalcMagnetization(int L, int (*lat)[L]); // calculate magnetization
double CalcEnergyOverall(int L, int (*lat)[L]); // calculate energy for overall lattice
double CalcEnergySite(int i, int j, int L, int (*lat)[L]); // calculate energy for particular site

void MonteCarlo(int L, int N, double T, Solution *s); // Monte Carlo simulation

#endif
