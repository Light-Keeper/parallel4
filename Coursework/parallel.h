#include <mpi.h>

class Matrix;
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

	static ParallelMatrixMultiplication *Instance()
	{
		static ParallelMatrixMultiplication *instance = new ParallelMatrixMultiplication();
		return instance;
	}
};