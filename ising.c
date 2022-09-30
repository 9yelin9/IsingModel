// ising.c : calculate 2D ising model

#include "ising.h"

FILE* OpenFile(char *fs, char *ftype) {
	FILE *f;

	if((f = fopen(fs, ftype)) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1); }

	return f;
}

void GenIdx(int N, int *idx) {
	int i, a, b, tmp;
	
	for(i=0; i<N; i++) idx[i] = i;

	for(i=0; i<100; i++) {
		a = (int)(10 * ran2(&seed)) % N;
		b = (int)(10 * ran2(&seed)) % N;

		if(a != b) {
			tmp    = idx[a];
			idx[a] = idx[b];
			idx[b] = tmp;
		}
	}
}

void GenName(int L, char *data_type, char *fs) {
	sprintf(fs,\
			"output/%s_MCS%d_L%d.txt",\
			data_type, (int)log10(MCS), L);
}

void InitLatRand(int L, int (*lat)[L]) {
	int i, j;

	for(i=0; i<L; i++) {
		for(j=0; j<L; j++) {
			if(ran2(&seed) > 0.5) lat[i][j] = 1;
			else lat[i][j] = -1;
		}
	}
}

void InitLatBin(int dec, int L, int (*lat)[L]) {
	int i, j;	

	for(i=0; i<L; i++) {
		for(j=0; j<L; j++) {
			lat[(L-1) - i][(L-1) - j] = 2 * (dec % 2) - 1;
			dec /= 2;
		}
	}
}

void InitMonteCarlo(int L, int N, double T, int (*lat)[L]) {
	int itr, i, I, J, idx[N];
	double e0, e1;

	for(itr=0; itr<EQS; itr++) {
		GenIdx(N, idx);

		for(i=0; i<N; i++) {
			I = idx[i] / L;
			J = idx[i] % L;

			e0 = CalcEnergySite(I, J, L, lat);
			lat[I][J] *= -1;
			e1 = CalcEnergySite(I, J, L, lat);
							
			if(e1 - e0 > 0 && ran2(&seed) > exp((e0 - e1) / T)) lat[I][J] *= -1;
		}
	}
}

double CalcMagnetization(int L, int (*lat)[L]) {
	int i, j;
	double m = 0;

	for(i=0; i<L; i++) {
		for(j=0; j<L; j++) {
			m += lat[i][j];
		}
	}

	return m;
}

double CalcEnergyOverall(int L, int (*lat)[L]) {
	int i, j, up, left;
	double e = 0;

	for(i=0; i<L; i++) {
		for(j=0; j<L; j++) {
			up    = i > 0 ? lat[i-1][j] : lat[L-1][j];
			left  = j > 0 ? lat[i][j-1] : lat[i][L-1];
			
			e += lat[i][j] * (up + left);
		}
	}

	return -e;
}

double CalcEnergySite(int i, int j, int L, int (*lat)[L]) {
	int up, dn, left, right;
	double e;

	up     = i > 0   ? lat[i-1][j] : lat[L-1][j];
	dn     = i < L-1 ? lat[i+1][j] : lat[0][j];
	left   = j > 0   ? lat[i][j-1] : lat[i][L-1];
	right  = j < L-1 ? lat[i][j+1] : lat[i][0];

	e = lat[i][j] * (up + dn + left + right);

	return -e;
}

void MonteCarlo(int L, int N, double T, MCObservable *mco) {
	int itr, i, I, J, lat[L][L], idx[N];
	double m, e, e0, e1;

	InitLatRand(L, lat);
	InitMonteCarlo(L, N, T, lat);
	
	mco->m = mco->ma = mco->m2 = mco->m4 = mco->e = mco->e2 = 0;
	
	for(itr=0; itr<MCS; itr++) {
		GenIdx(N, idx);

		for(i=0; i<N; i++) {
			I = idx[i] / L;
			J = idx[i] % L;

			e0 = CalcEnergySite(I, J, L, lat);
			lat[I][J] *= -1;
			e1 = CalcEnergySite(I, J, L, lat);
							
			if(e1 - e0 > 0 && ran2(&seed) > exp((e0 - e1) / T)) lat[I][J] *= -1;
		}

		m = CalcMagnetization(L, lat);
		e = CalcEnergyOverall(L, lat);

		mco->m  += m;
		mco->ma += abs(m);
		mco->m2 += pow(m, 2);
		mco->m4 += pow(m, 4);
		mco->e  += e;
		mco->e2 += pow(e, 2);
	}

	mco->m  /= MCS;
	mco->ma /= MCS;
	mco->m2 /= MCS;
	mco->m4 /= MCS;
	mco->e  /= MCS;
	mco->e2 /= MCS;
}

void Test() {
	int L = 2, N = L*L;

	int itr, i, j, lat[L][L];
	double e;

	for(itr=0; itr<pow(2, N); itr++) {
		InitLatBin(itr, L, lat);

		e = CalcEnergyOverall(L, lat);

		for(i=0; i<L; i++) {
			for(j=0; j<L; j++) {
				printf("%2d\t", lat[i][j]);
			}
		}
		printf("%f\n", e);
	}	
}

int main(int argc, char *argv[]) {
	if(argc != 2) {
		printf("%s <L> : make L*L Ising model\n", argv[0]);
		exit(1);
	}

	int L = atoi(argv[1]), N = L*L, t;
	double T;
	MCObservable mco;

	FILE *fe, *fc, *fm, *fx, *fu;
	char fes[64], fcs[64], fms[64], fxs[64], fus[64];

	GenName(L, "energy", fes);
	GenName(L, "heat_capacity", fcs);
	GenName(L, "abs_magnetization", fms);
	GenName(L, "mag_susceptibility", fxs);
	GenName(L, "cumulant", fus);

	fe = OpenFile(fes, "w");
	fc = OpenFile(fcs, "w");
	fm = OpenFile(fms, "w");
	fx = OpenFile(fxs, "w");
	fu = OpenFile(fus, "w");

	for(t=1; t<T_MAX; t++) {
		T = 0.1 * (T_MAX - t);

		MonteCarlo(L, N, T, &mco);

		fprintf(fe, "%f\t%f\n", T, ENERGY / N);
		fprintf(fc, "%f\t%f\n", T, HEAT_CAPACITY(T) / N);
		fprintf(fm, "%f\t%f\n", T, ABS_MAGNETIZATION / N);
		fprintf(fx, "%f\t%f\n", T, MAG_SUSCEPTIBILITY(T) / N);
	}

	for(t=200; t<300; t++) {
		T = 0.01 * (300 - t);

		MonteCarlo(L, N, T, &mco);

		fprintf(fu, "%f\t%f\n", T, CUMULANT / N);
	}

	fclose(fe);
	fclose(fc);
	fclose(fm);
	fclose(fx);
	fclose(fu);

	return 0;
}
