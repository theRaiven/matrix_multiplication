#include "Multiply.h"

int main()
{
	setlocale(LC_ALL, "rus");
	int n = 2500;
	int threads[]{ 2, 16 };
	double* A = new double[n * n];
	double* B = new double[n * n];
	double* C = new double[n * n];

	FillVector(A, n * n);
	FillVector(B, n * n);


	cout << "===== Последовательная версия =====" << endl;

	double start_seq = omp_get_wtime();
	MultiplicationMatrixSeq(A, B, C, n);
	double end_seq = omp_get_wtime();
	cout << endl << "Время выполнения: " << (end_seq - start_seq)
		<< " секунд" << endl;
	cout << "Норма матрицы C: " << NormMatrix(C, n) << endl;

	//============================================================
	for (int j = 0; j < 2; j++)
	{

		cout << "\n===== Параллельная версия ("
			<< threads[j] << " потока) =====" << endl;

		double start_par = omp_get_wtime();
		MultiplicationMatrixPar(A, B, C, n, threads[j]);
		double end_par = omp_get_wtime();
		cout << endl << "Время выполнения: " <<
			(end_par - start_par) << " секунд" << endl;
		cout << "Норма матрицы C: " << NormMatrix(C, n) << "\tУскорение: " << (end_seq - start_seq) / (end_par - start_par) << endl;
	}
	delete[] A;
	delete[] B;
	delete[] C;
}
