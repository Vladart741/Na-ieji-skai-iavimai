#include "pch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <chrono>
#include <mpi.h>

using namespace std;

int numDP = 1000;      // Vietoviu skaicius (demand points, max 10000)
int numPF = 10;        // Esanciu objektu skaicius (preexisting facilities)
int numCL = 22;        // Kandidatu naujiems objektams skaicius (candidate locations)
int numX = 2;         // Nauju objektu skaicius

double **demandPoints; // Geografiniai duomenys


//=============================================================================

void loadDemandPoints();
double HaversineDistance(double* a, double* b);
double evaluateSolution(int *X);
int increaseX(int *X, int index, int maxindex);

//=============================================================================

int main(int argc, char *argv[]) {

	int id, numProcs;
	MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &id);
	MPI_Comm_size(comm, &numProcs);

	if (id == 0) {
		cout << endl << "Algoritmas naudojant " << numProcs << " Threads:" << endl << endl;
	}
start:
	auto start1 = std::chrono::high_resolution_clock::now();
	numX++;
	if (id == 0) {
		cout << "Kai numX = " << numX << endl;
	}

	


	loadDemandPoints();             // Nuskaitomi duomenys
	MPI_Barrier(comm);
	// Sudarom pradini sprendini: [0, 1, 2, 3, ...]
	int *X = new int[numX];
	int *bestX = new int[numX];
	for (int i = 0; i < numX; i++) {
		X[i] = i;
		bestX[i] = i;
	}
	double u = evaluateSolution(X);
	double bestU = u;
	int r;
	//----- Pagrindinis ciklas ------------------------------------------------


	while (true) {
		if (increaseX(X, numX - 1, numCL)) {
			u = evaluateSolution(X);
			if (u > bestU) {
				bestU = u;
				for (int i = 0; i < numX; i++) bestX[i] = X[i];
			}
		}
		else break;
	}
	//----- Rezultatu spausdinimas --------------------------------------------
	MPI_Barrier(comm);
	if (id == 0) {

		auto finish1 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed1 = finish1 - start1;

		cout << "Geriausias sprendinys: ";
		for (int i = 0; i < numX; i++) cout << bestX[i] << " ";
		cout << "(" << bestU << ")" << endl;
		cout << "Skaiciavimo trukme: " << elapsed1.count() << endl;
	}
	MPI_Barrier(comm);
	if (numX < 5) goto start;
	MPI_Finalize();
	return 0;
}

//=============================================================================

void loadDemandPoints() {
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double*[numDP];
	for (int i = 0; i < numDP; i++) {
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

double HaversineDistance(double* a, double* b) {
	double dlon = fabs(a[0] - b[0]);
	double dlat = fabs(a[1] - b[1]);
	double aa = pow((sin((double)dlon / (double)2 * 0.01745)), 2) + cos(a[0] * 0.01745) * cos(b[0] * 0.01745) * pow((sin((double)dlat / (double)2 * 0.01745)), 2);
	double c = 2 * atan2(sqrt(aa), sqrt(1 - aa));
	double d = 6371 * c;
	return d;
}


//=============================================================================

double evaluateSolution(int *X) {
	int id, numProcs;

	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &id);
	MPI_Comm_size(comm, &numProcs);

	int defaultNumTasks = numDP / numProcs;
	int numtasks;
	int remainder = numDP % numProcs;
	int start;

	if (remainder == 0) {
		numtasks = numDP / numProcs;
		start = numtasks * id;
	}
	else {
		if (id < remainder) {
			numtasks = defaultNumTasks + 1;
			start = numtasks * id;
		}
		else {
			numtasks = defaultNumTasks;
			start = numtasks * id + remainder;
		}
	}

	MPI_Barrier(comm);

	double U = 0;
	int bestPF;
	int bestX;
	double d;
	for (int i = start; i < (start + numtasks); i++) {
		bestPF = 1e5;
		for (int j = 0; j < numPF; j++) {
			d = HaversineDistance(demandPoints[i], demandPoints[j]);
			if (d < bestPF) bestPF = d;
		}
		bestX = 1e5;
		for (int j = 0; j < numX; j++) {
			d = HaversineDistance(demandPoints[i], demandPoints[X[j]]);
			if (d < bestX) bestX = d;
		}
		if (bestX < bestPF) U += demandPoints[i][2];
		else if (bestX == bestPF) U += 0.3*demandPoints[i][2];
	}
	MPI_Barrier(comm);

	double *rbuf;
	if (id == 0) {
		rbuf = (double *)malloc(numProcs * sizeof(double));
	}
	else {
		rbuf = (double *)malloc(0);
	}

	MPI_Barrier(comm);
	MPI_Gather(&U, 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, comm);

	U = 0;

	if (id == 0) {
		for (int i = 0; i < numProcs; i++) {
			U = U + rbuf[i];
		}
	}

	MPI_Barrier(comm);

	return U;
}

//=============================================================================

int increaseX(int *X, int index, int maxindex) {
	if (X[index] + 1 < maxindex - (numX - index - 1)) {
		X[index]++;
	}
	else {
		if ((index == 0) && (X[index] + 1 == maxindex - (numX - index - 1))) {
			return 0;
		}
		else {
			if (increaseX(X, index - 1, maxindex)) X[index] = X[index - 1] + 1;
			else return 0;
		}
	}
	return 1;
}
