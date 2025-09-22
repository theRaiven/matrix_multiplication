#include "Multiply.h"

void FillVector(double*& vect, int n)
{
	srand(static_cast<unsigned>(time(0)));

	for (int i = 0; i < n; i++)
	{
		vect[i] = rand() % 10;
	}
}
double MultiplicationVectorSeq(double* x, double* y, int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		sum += x[i] * y[i];
	}
	return sum;
}

double MultiplicationVectorPar(double* x, double* y, int n, int threads)
{
	double sum = 0.0;
	omp_set_num_threads(threads);

#pragma omp parallel
	{
		// один поток выведет реальное количество
		if (omp_get_thread_num() == 0)
		{
			std::cout << "Запущено потоков: "
				<< omp_get_num_threads() << std::endl;
		}

#pragma omp for reduction(+:sum)
		for (int i = 0; i < n; i++)
		{
			sum += x[i] * y[i];
		}
	}

	return sum;
}

void MultiplicationMatrixSeq(double* A, double* B, double* C, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double sum = 0.0;
			for (int k = 0; k < n; k++)
			{
				sum += A[i * n + k] * B[k * n + j];
			}
			C[i * n + j] = sum;
		}
	}
}
void MultiplicationMatrixPar(double* A, double* B, double* C, int n, int threads)
{
	omp_set_num_threads(threads);

#pragma omp parallel
	{

		if (omp_get_thread_num() == 0)
		{
			std::cout << "Запущено потоков: "
				<< omp_get_num_threads() << std::endl;
		}

#pragma omp for collapse(2) schedule(dynamic, 16)
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				double sum = 0.0;
				for (int k = 0; k < n; k++)
				{
					sum += A[i * n + k] * B[k * n + j];
				}
				C[i * n + j] = sum;
			}
		}
	}
}
double NormMatrix(double* C, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n * n; i++)
	{
		sum += C[i] * C[i];
	}
	return sqrt(sum);
}

void FillUpperTrianglar(double* A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j < i)
			{
				A[i * n + j] = 0.0;
			}
			else if (i == j)
			{
				A[i * n + j] = n + i;
			}
			else
			{
				A[i * n + j] = rand() % 10 + 1;
			}
		}
	}
}

void MatVec(double* A, double* x, double* b, int n)
{
	for (int i = 0; i < n; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < n; j++)
		{
			sum += A[i * n + j] * x[j];
		}
		b[i] = sum;
	}
}

void BackSubstitutionSeq(double* A, double* b, double* x, int n)
{
	for (int i = n - 1; i >= 0; i--)
	{
		double sum = b[i];
		for (int j = i + 1; j < n; j++)
		{
			sum -= A[i * n + j] * x[j];
		}
		x[i] = sum / A[i * n + i];
	}
}

void BackSubstitutionPar(double* A, double* b, double* x, int n, int threads)
{
	omp_set_num_threads(threads);

    for (int i = n - 1; i >= 0; i--)
    {
        double sum = b[i];

        #pragma omp parallel for reduction(-:sum)
        for (int j = i + 1; j < n; j++)
        {
            sum -= A[i * n + j] * x[j];
        }

        x[i] = sum / A[i * n + i];
    }
}
double Norm(double* y, int n)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
	{
		sum += y[i] * y[i];
	}
	return sqrt(sum);
}
void MatVecSeq(double* A, double* x, double* y, int n)
{
	for (int i = 0; i < n; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < n; j++)
		{
			sum += A[i * n + j] * x[j];
		}
		y[i] = sum;
	}
}

void MatVecPar(double* A, double* x, double* y, int n, int threads)
{
#pragma omp parallel num_threads(threads)
	{
		if (omp_get_thread_num() == 0)
			std::cout << "Запущено потоков: " << omp_get_num_threads() << "\t";

#pragma omp for
		for (int i = 0; i < n; i++)
		{
			double sum = 0.0;
			for (int j = 0; j < n; j++)
			{
				sum += A[i * n + j] * x[j];
			}
			y[i] = sum;
		}
	}
}

void FillMatrixVector(double* A, double* x, int n)
{
	for (int i = 0; i < n; i++)
	{
		{
			x[i] = rand() % 10 + 1;
		}
		for (int j = 0; j < n; j++)
		{
			A[i * n + j] = rand() % 10 + 1;
		}
	}
}

