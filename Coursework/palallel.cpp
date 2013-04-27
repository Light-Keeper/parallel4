#pragma omp

#include "stdio.h"
#include "mpi.h"
#include "../Matrix/Matrix.h"

class ParallelMatrixMultiplication
{
	static const int EVENT_MUL = 0;
	static const int EVENT_EXIT = 1;
	static const int TAG_CMD = 3;
	static const int TAG_DATA_1 = 4;
	static const int TAG_DATA_2 = 5;
	static const int TAG_DATA_3 = 6;

	int CurrentNode;
	int NumberOfNodes;
	MPI_Status status;

	double *a;
	double *b;
	double *c;

	void DispatchEvents();
	void MulHelper();
public:
	bool Init(int argc, char **argv);
	bool Finalize();
	Matrix Mul(const Matrix &x, const Matrix &y);
};


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
	for (int i=1; i <= NumberOfNodes; i++){
		MPI_Recv(&info, sizeof(info), MPI_BYTE, i, TAG_DATA_1, MPI_COMM_WORLD, &status);
		MPI_Recv(&result[info.offset][0], info.rows* info.A_width, MPI_DOUBLE, i, TAG_DATA_2, MPI_COMM_WORLD, &status);
	}

	return result;
}

main(int argc, char **argv){
	struct timeval start, stop;
	int numberOfTasks,
		mtype,
		taskID,
		numberOfWorkers,
		source,
		destination,
		rows,
		averageRow,
		extra,
		offset,i,j,k;
	//first initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskID);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfTasks);

	numberOfWorkers = numberOfTasks-1;

	//---------------------------- master ----------------------------//
	if (taskID == 0) {
		for (i=0; i<N; i++) {
			for (j=0; j<N; j++) {   
				a[i][j]= 1.0;
				b[i][j]= 2.0;
			}
		}

		/* send matrix data to the worker tasks */

		gettimeofday(&start, 0);

		averageRow = N/numberOfWorkers; //average rows per worker
		extra = N % numberOfWorkers;  //extra rows
		offset = 0;

		for (destination=1; destination<=numberOfWorkers; destination++)
		{       

			if(destination <= extra){
				rows = averageRow+1;
			}else{
				rows = averageRow;
			}
			mtype = 1;
			MPI_Send(&offset, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, destination, mtype, MPI_COMM_WORLD);
			MPI_Send(&a[offset][0], rows*N, MPI_DOUBLE,destination,mtype, MPI_COMM_WORLD);
			MPI_Send(&b, N*N, MPI_DOUBLE, destination, mtype, MPI_COMM_WORLD);
			offset = offset + rows;
		}

		/* wait for results from all worker tasks */
		for (i=1; i<=numberOfWorkers; i++){
			mtype = 2;
			source = i;
			MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&c[offset][0], rows*N, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
		}
		gettimeofday(&stop, 0);


		printf("Upper Left = %6.2f    Upper Right = %6.2f \n",c[0][0],c[0][N-1]);
		printf("Lower Left = %6.2f    Lower Right = %6.2f \n",c[N-1][0],c[N-1][N-1]);

		fprintf(stdout,"Time = %.6f\n\n",(stop.tv_sec+stop.tv_usec*1e-6)-(start.tv_sec+start.tv_usec*1e-6));


	} 

	/*---------------------------- worker----------------------------*/
	if (taskID > 0) {
		source = 0;
		mtype = 1;
		MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&a, rows*N, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&b, N*N, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);

		/* Matrix multiplication */
		for (k=0; k<N; k++)
			for (i=0; i<rows; i++) {
				c[i][k] = 0.0;
				for (j=0; j<N; j++)
					c[i][k] = c[i][k] + a[i][j] * b[j][k];
			}

			mtype = 2;
			MPI_Send(&offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
			MPI_Send(&c, rows*N, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
	}  

	MPI_Finalize();
} 