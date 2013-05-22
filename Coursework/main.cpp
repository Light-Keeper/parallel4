#include <cstdio>
#include <algorithm>
#include <assert.h>
#include "Matrix.h"
#include "parallel.h"
#include <cmath>

#define MAX_ERROR 0.000001
#define STEP 0.1

// callbacks
void set0(double *x)
{
	*x = 0;
}

void set_X_AbsX(double *x)
{
	* x =  *x * fabs(*x);
}

void SetRandomNumber(double *x)
{
	*x = rand() / 3;
}



int MainThread()
{
	FILE *f = fopen("t.txt", "r");
	if (f == NULL)
	{
		printf("input file not found\n");
		return 0;
	}

	double runtime = - MPI_Wtime();
	
	// test big matrix multiplication 
	if ( 0 )
	{
		Matrix m1(1000, 1000);
		Matrix m2(1000, 400);
		m1.foreach( SetRandomNumber );
		m2.foreach( SetRandomNumber );

		Matrix m3  = m1 * m2;
		runtime += MPI_Wtime();
		printf("Calculation finished in %lf sec", runtime);
		return 0;
	}


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

	runtime += MPI_Wtime();

	printf("Calculation finished in %lf sec", runtime);
	fclose(f);	
	return 0;
}



int main(int argc, char *argv[])
{
	// initialize MPI
	ParallelMatrixMultiplication::Instance()->Init(argc, argv);
	MainThread();
	ParallelMatrixMultiplication::Instance()->Finalize();
	return 0;
}

// this method initializes MPI. It returns only for process with rank = 0.
// for other it starts DispatchEvents.
bool ParallelMatrixMultiplication::Init(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &CurrentNode);
	MPI_Comm_size(MPI_COMM_WORLD, &NumberOfNodes);
	
	// suppouse we have matrizes less then 1000 * 1000
	a = (double *)malloc( sizeof(double) * 1000 * 1000 );
	b = (double *)malloc( sizeof(double) * 1000 * 1000 );
	c = (double *)malloc( sizeof(double) * 1000 * 1000 );
	assert(a != NULL);
	assert(b != NULL);
	assert(c != NULL);

	if (CurrentNode == 0) return true; 
	DispatchEvents();
	return true;
}

// finalize MPI, free memory and exit
bool ParallelMatrixMultiplication::Finalize()
{
	free( a );
	free( b );
	free( c );
		
	if (CurrentNode == 0)
	{
		// send EVENT_EXIT to all other processes
		for (int i = 1; i < NumberOfNodes; i++)
		{
			int code = EVENT_EXIT;
			MPI_Send(&code, 1, MPI_INT, i, TAG_CMD, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize();
	return true;
}

//waits for commands and starts appropriate routine
void ParallelMatrixMultiplication::DispatchEvents()
{
	while ( 1 )
	{
		// get command from process with rank = 0
		int cmd = 0;
		MPI_Recv(&cmd, 1, MPI_INT, 0, TAG_CMD, MPI_COMM_WORLD, &status);
		switch (cmd)
		{
		case EVENT_MUL:
			// process 0 wants us to help him to multiply matrices.
			MulHelper();
			break;
		case EVENT_EXIT:
			// all is done
			Finalize();
			exit( 0 );
		default:
			break;
		}
	}
}

//receive part of first matrix, entire second, and multiply it.
void ParallelMatrixMultiplication::MulHelper()
{
	int source = 0;
	MatrixInfo info;

	// recceive
	MPI_Recv(&info,	sizeof(info), MPI_BYTE	, source, TAG_DATA_1, MPI_COMM_WORLD, &status);
	MPI_Recv(a,   info.rows * info.A_width, MPI_DOUBLE	, source, TAG_DATA_2, MPI_COMM_WORLD, &status);
	MPI_Recv(b,   info.B_height * info.B_width, MPI_DOUBLE	, source, TAG_DATA_3, MPI_COMM_WORLD, &status);

	// multiply
	for ( int i = 0; i < info.rows; i++)
		for( int j = 0; j < info.B_width; j++)
		{
			double result;
			double * _a = &(a[i * info.A_width]) - 1;
			double * _b = &(b[j]);
			result = 0.0;
			
			for ( int t = 0; t < info.A_width; t++)
			{
				result += *(++_a) * *_b;
				_b += info.B_width; // next row
			}

			c[i * info.B_width + j] = result;
		}

	// send result back
	MPI_Send(&info, sizeof(info), MPI_BYTE, 0, TAG_DATA_1, MPI_COMM_WORLD);
	MPI_Send(c, info.rows * info.B_width, MPI_DOUBLE, 0, TAG_DATA_2, MPI_COMM_WORLD);
}

// multiply 2 matrices. call it from main thread. 
Matrix ParallelMatrixMultiplication::Mul(const Matrix &x, const Matrix &y)
{
	Matrix result(x.n, y.m);
	
	if (NumberOfNodes == 1)
	{
		// only 1 node, make it just here.
		for (int i = 0; i < x.n; i++) 
			for(int j = 0; j < y.m; j++)
				{
					result[i][j] = 0;
					for (int t = 0; t < x.m; t++) 
						result[i][j] += x[i][t] * y[t][j];
				}
			return result;
	}

	int averageRow = x.n / (NumberOfNodes - 1); //average rows per worker
	int	extra = x.n % (NumberOfNodes - 1);  //extra rows
	MatrixInfo info;
	info.A_height = x.n;
	info.A_width = x.m;
	info.B_height = y.n;
	info.B_width = y.m;
	info.offset = 0;
	
	// send to each node part of job
	for (int destination=1; destination < NumberOfNodes; destination++)
	{       
		info.rows = averageRow  + (destination <= extra);
		int cmd = EVENT_MUL;
		MPI_Send(&cmd, 1, MPI_INT, destination, TAG_CMD, MPI_COMM_WORLD);	
		MPI_Send(&info, sizeof(info), MPI_BYTE, destination, TAG_DATA_1, MPI_COMM_WORLD);
		MPI_Send(&(x[info.offset][0]), info.rows*x.m, MPI_DOUBLE,destination, TAG_DATA_2, MPI_COMM_WORLD);
		MPI_Send(&(y[0][0]), y.n * y.m, MPI_DOUBLE, destination, TAG_DATA_3, MPI_COMM_WORLD);
		info.offset = info.offset + info.rows;
	}

	// wait for results from all worker tasks 
	for (int i=1; i < NumberOfNodes; i++)
	{
		MPI_Recv(&info, sizeof(info), MPI_BYTE, i, TAG_DATA_1, MPI_COMM_WORLD, &status);
		MPI_Recv(&(result[info.offset][0]), info.rows* info.B_width, MPI_DOUBLE, i, TAG_DATA_2, MPI_COMM_WORLD, &status);
	}
	return result;
}
