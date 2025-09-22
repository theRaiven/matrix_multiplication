#pragma once
#include <iostream>
#include <fstream>
#include <cstdlib>   // для rand() и srand()
#include <ctime>     // для time()
#include <omp.h>
#include <iomanip>
using namespace std;

double MultiplicationVectorSeq(double* x, double* y, int n);
double MultiplicationVectorPar(double* x, double* y, int n, int threads);

void MultiplicationMatrixSeq(double* A, double* B, double* C, int n);
void MultiplicationMatrixPar(double* A, double* B, double* C, int n, int threads);
double NormMatrix(double* C, int n);
void FillVector(double*& vect, int n);

void BackSubstitutionPar(double* A, double* b, double* x, int n, int threads);
void BackSubstitutionSeq(double* A, double* b, double* x, int n);
void MatVec(double* A, double* x, double* b, int n);
void FillUpperTrianglar(double* A, int n);
double Norm(double* x, int n);

void FillMatrixVector(double* A, double* x, int n);
void MatVecSeq(double* A, double* x, double* y, int n);
void MatVecPar(double* A, double* x, double* y, int n, int threads);
