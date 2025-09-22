#include "Multiply.h"
void FirstEx();
void SecondEx();
void TrithEx();
void FourthEx();

int main()
{
	setlocale(LC_ALL, "rus");
	FirstEx();
}
void FirstEx()
{
	int n[]{ 10, 1000, 100000, 1000000 };
	int threads[]{ 2, 16 };
	for (int i = 0; i < 4; i++)
	{
		cout << "******** Тест для размерности " << n[i] << " ********" << endl;
		double* x = new double[n[i]];
		double* y = new double[n[i]];

		FillVector(x, n[i]);
		FillVector(y, n[i]);


		cout << "===== Последовательная версия =====" << endl;

		double start_seq{ omp_get_wtime() };

	double sumSeq{ MultiplicationVectorSeq(x, y, n[i]) };
	double end_seq{ omp_get_wtime() };
	cout << endl << "Время выполнения: "
		<< fixed << setprecision(7) << (end_seq - start_seq) <<
		" секунд" << endl;
	cout << "Результат: " << sumSeq << endl;


	//============================================================
	for (int j = 0; j < 2; j++)
	{
		cout << "\n===== Параллельная версия (" << threads[j] <<
			" потока) =====" << endl;
		double start_par = omp_get_wtime();

		double sumPar = MultiplicationVectorPar(x, y, n[i], threads[j]);
		double end_par = omp_get_wtime();
		cout << endl << "Время выполнения: "
			<< fixed << setprecision(7) << (end_par - start_par) <<
			" секунд" << endl;
		cout << "Результат: " << sumPar << "\tУскорение: "
			<< (end_seq - start_seq) / (end_par - start_par) << endl;
	}


	delete[] x;
	delete[] y;
}
}

void SecondEx()
{
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
void TrithEx()
{
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
void FourthEx()
{
	int n = 2500;
	int threads[]{ 2, 16 };

	double* A = new double[n * n];
	double* x_true = new double[n];
	double* b = new double[n];
	double* x_seq = new double[n];
	double* x_par = new double[n];

	FillUpperTrianglar(A, n);
	FillVector(x_true, n);

	// b = A * x_true
	MatVec(A, x_true, b, n);

	cout << "===== Последовательная версия =====" << endl;
	double start_seq = omp_get_wtime();
	BackSubstitutionSeq(A, b, x_seq, n);
	double end_seq = omp_get_wtime();
	cout << fixed << setprecision(6);
	cout << endl << "Время выполнения: " << (end_seq - start_seq) << " секунд" << endl;

	//============================================================
	for (int j = 0; j < 2; j++)
	{
		cout << "\n===== Параллельная версия (" << threads[j] << " потока) =====" << endl;
		double start_par = omp_get_wtime();
		BackSubstitutionPar(A, b, x_par, n, threads[j]);
		double end_par = omp_get_wtime();

		cout << fixed << setprecision(6);
		cout << endl << "Время выполнения: " <<
			(end_par - start_par) << " секунд" << endl;
		cout << "Норма x_par  = " << Norm(x_par, n) << "\tУскорение: " << (end_seq - start_seq) / (end_par - start_par) << endl;
	}

	//============================================================
	ofstream fout("results.txt");
	fout << fixed << setprecision(6);
	fout << "x_true vs x_seq vs x_par\n";
	for (int i = 0; i < n; i++)
	{
		fout << x_true[i] << " " << x_seq[i] << " " << x_par[i] << "\n";
	}
	fout.close();

	cout << "\nНорма x_true = " << Norm(x_true, n) << endl;
	cout << "Норма x_seq  = " << Norm(x_seq, n) << endl;
	cout << "Норма x_par  = " << Norm(x_par, n) << endl;

	delete[] A;
	delete[] x_true;
	delete[] b;
	delete[] x_seq;
	delete[] x_par;

}