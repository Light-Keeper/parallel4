#include <cstdio>
#include <algorithm>
#include <assert.h>
#include "Matrix.h"

#define MAX_ERROR 0.00001
#define STEP 0.001

void set0(double *x)
{
	*x = 0;
}

void set_X_AbsX(double *x)
{
	* x = *x * fabs(*x);
}

int main(int argc, char *argv[])
{
	FILE *f = fopen("t.txt", "r");
	ASSERT(f != NULL);
	
	Matrix Ax, Ay, Sx, Sy, R, Kx, Ky, H;
	
	Ax = Matrix::ReadFromFile(f);
	Ay = Matrix::ReadFromFile(f);
	
	Sx = Matrix::ReadFromFile(f);
	Sy = Matrix::ReadFromFile(f);
	
	Kx = Matrix::ReadFromFile(f);
	Ky = Matrix::ReadFromFile(f);
	
	R = Matrix::ReadFromFile(f);
	H = Matrix::ReadFromFile(f);
		
	fclose(f);
	
	Matrix W  = Ax.Inverse() * Ay;	
	Matrix TP = (Sy * Ky - Sx * Kx * W).Inverse() * (Sx & Sy);	
	Matrix Y(TP.n, 1);
	Matrix X(W.n,  1);
	
	Y.foreach( set0 );
	X.foreach( set0 );

	Matrix LastZ = (X.Transpose() & Y.Transpose()).Transpose();
	Matrix Z(LastZ);

	do 
	{
		LastZ = Z;

		Y = Y + TP * STEP * (H - R * Z);
		X = W * Y * -1;

		Z = (X.Transpose() & Y.Transpose()).Transpose();
		Z.foreach( set_X_AbsX );

	} while ( Z.SqrDifference( LastZ ) > sqr( MAX_ERROR ) );

	f = fopen("out2.txt", "w");
	ASSERT( f );
	
	fprintf(f, "W(%d, %d) = \n", W.n, W.m);
	W.PrintToFile( f );
	fprintf(f, "TP(%d, %d) = \n", TP.n, TP.m);
	TP.PrintToFile( f );

	for (int i = 0; i < X.n; i++)
		fprintf(f, "X[ %d ] = %lf\n", i, X[i][0]);
	
	for (int i = 0; i < Y.n; i++)
		fprintf(f, "Y[ %d ] = %lf\n", i, Y[i][0]);

	fclose(f);	
	return 0;
}




bool ParallelMatrixMultiplication::Init(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &CurrentNode);
	MPI_Comm_size(MPI_COMM_WORLD, &NumberOfNodes);
	
	a = (double *)malloc( sizeof(double) * 500 * 500 );
	b = (double *)malloc( sizeof(double) * 500 * 500 );
	c = (double *)malloc( sizeof(double) * 500 * 500 );

	if (CurrentNode == 0) return true;
	DispatchEvents();
	return true;
}

bool ParallelMatrixMultiplication::Finalize()
{
	MPI_Finalize();
	free( a );
	free( b );
	free( c );

	return true;
}


void ParallelMatrixMultiplication::DispatchEvents()
{
	while ( 1 )
	{
		int cmd = 0;
		MPI_Recv(&cmd, 1, MPI_INT, 0, TAG_CMD, MPI_COMM_WORLD, &status);
		switch (cmd)
		{
		case EVENT_MUL:
			MulHelper();
			break;
		case EVENT_EXIT:
			Finalize();
			return;
		default:
			break;
		}
	}
}

struct MatrixInfo
{
	int A_width;
	int A_height;
	int B_width;
	int B_height;
	int offset;
	int rows;
};

void ParallelMatrixMultiplication::MulHelper()
{
	int source = 0;
	MatrixInfo info;

	MPI_Recv(&info,	sizeof(info), MPI_BYTE	, source, TAG_DATA_1, MPI_COMM_WORLD, &status);
	MPI_Recv(&a,   info.rows * info.A_width, MPI_DOUBLE	, source, TAG_DATA_2, MPI_COMM_WORLD, &status);
	MPI_Recv(&b,   info.B_height * info.B_width, MPI_DOUBLE	, source, TAG_DATA_3, MPI_COMM_WORLD, &status);
			
	for (int i = 0; i < info.rows; i++)
		for(int j = 0; j < info.B_width; j++)
		{
			double *result = &c[i * info.B_width + j];
			double * _a = &a[i * info.A_width];
			double * _b = &b[j];
			*result = 0.0;
			
			for (int t = 0; t < info.A_width; t++)
			{
				*result += *_a * *_b;
				_a++; // следующий столбец
				_b += info.B_width; // следующая строка
			}
		}

	MPI_Send(&info, sizeof(info), MPI_INT, 0, TAG_DATA_1, MPI_COMM_WORLD);
	MPI_Send(&c, info.rows * info.A_width, MPI_DOUBLE, 0, TAG_DATA_2, MPI_COMM_WORLD);
}

Matrix ParallelMatrixMultiplication::Mul(const Matrix &x, const Matrix &y)
{
	Matrix result(x.n, y.m);

	int averageRow = x.n / (NumberOfNodes - 1); //average rows per worker
	int	extra = x.n % (NumberOfNodes - 1);  //extra rows
	MatrixInfo info;
	info.A_height = x.n;
	info.A_width = x.m;
	info.B_height = x.n;
	info.B_width = x.m;
	info.offset = 0;
	
	for (int destination=1; destination <= NumberOfNodes; destination++)
	{       
		info.rows = averageRow  + (destination <= extra);
		MPI_Send(EVENT_MUL, 1, MPI_INT, destination, TAG_CMD, MPI_COMM_WORLD);	
		MPI_Send(&info, sizeof(info), MPI_BYTE, destination, TAG_DATA_1, MPI_COMM_WORLD);
		MPI_Send(&x[info.offset][0], info.rows*x.m, MPI_DOUBLE,destination, TAG_DATA_2, MPI_COMM_WORLD);
		MPI_Send(&y[0][0], y.n * y.m, MPI_DOUBLE, destination, TAG_DATA_3, MPI_COMM_WORLD);
		info.offset = info.offset + info.rows;
	}

		/* wait for results from all worker tasks */
	for (int i=1; i <= NumberOfNodes; i++)
	{
		MPI_Recv(&info, sizeof(info), MPI_BYTE, i, TAG_DATA_1, MPI_COMM_WORLD, &status);
		MPI_Recv(&result[info.offset][0], info.rows* info.A_width, MPI_DOUBLE, i, TAG_DATA_2, MPI_COMM_WORLD, &status);
	}

	return result;
}
