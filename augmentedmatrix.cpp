#include "augmentedmatrix.h"

AugmentedMatrix::AugmentedMatrix(const Matrix& m, const Vector& v) : matrix(m), vector(v) {
}

AugmentedMatrix::~AugmentedMatrix() {

}

Vector AugmentedMatrix::cramersRule() const {
	double mainDeterminant = matrix.getDeterminant();
	if (mainDeterminant == 0)
		throw QString("Cannot solve this way. Determinant == 0");
	int row = matrix.rows();
	int col = matrix.cols();
	double** array;
	double* result = new double[row];
	for (int i = 0; i < row; i++) {
        array = matrix.getMatrix();
        Matrix::swapVectCol(array, vector, i);
		result[i] = matrix.getDeterminant(array, row, col) / mainDeterminant;
	}
	return Vector(result, row);
}

Vector AugmentedMatrix::inverseMatrixRule() const {
	double mainDeterminant = matrix.getDeterminant();
	if (mainDeterminant == 0)
		throw QString("Cannot solve this way");
	Matrix inversibleMatrix = this->matrix.inversibleMatrix();
	Vector vect = inversibleMatrix.multiply(vector);
	return vect;
}

Vector AugmentedMatrix::gaussianElimination() const {
	if (matrix.rows() != vector.size())
		throw QString("Cannot create matrix from this vector and matrix.");
	int row = matrix.rows();
	int col = matrix.cols() + 1;
	double** augMatrix = new double*[row];
	for (int i = 0; i < row; i++) {
		augMatrix[i] = new double[col];
		for (int j = 0; j < col - 1; j++) {
            augMatrix[i][j] = matrix(i, j);
		}
		augMatrix[i][col - 1] = vector.get(i);
	}
	this->matrix.forwardElimination(augMatrix, row, col);
	this->matrix.backwardElimination(augMatrix, row, col);
	double* result = new double[row];
	for (int i = 0; i < row; i++)
		result[i] = augMatrix[i][col - 1];
	return Vector(result, row);
}

Vector AugmentedMatrix::simpleIterationElimination(double eps) const {
	int row = rows();
	int col = cols() + 1;
	double** matrix = AugmentedMatrix::getAugmentedMatrix(this->matrix, this->vector);
	Matrix::normalizeMatrix(matrix, row, col);
	double a8 = Matrix::maximumAbsoluteRowSumNorm(matrix, row, col - 1);
	double a1 = Matrix::maximumAbsoluteColSumNorm(matrix, row, col - 1);
	double a2 = Matrix::spectralNorm(matrix, row, col - 1);
	double q = min(3, a8, a1, a2);
	if (q > 1)
		throw QString("Cannot find solution this way. (A q > 1)");
	double* result = new double[row];
	double* remember = new double[row];
	for (int i = 0; i < row; i++) {
		remember[i] = vector.get(i);
		result[i] = vector.get(i);
	}

	double d8 = Vector::maximumAbsoluteRowSumNorm(result, row);
	int N = (log(eps * (1 - q) / d8) / log(q)) + 1;
	for (int k = 0; k < N; k++) {
		for (int i = 0; i < row; i++) {
			remember[i] = result[i];
			result[i] = 0;
			for (int j = 0; j < col - 1; j++) {
				result[i] += matrix[i][j] * remember[j];
			}
			result[i] += matrix[i][col - 1];
		}
	}
	Matrix::freeMatrix(matrix, row);
	Vector::freeVector(remember);
	return Vector(result, row);
}

Vector AugmentedMatrix::seidelIterationElimination(double eps) const {
	int row = rows();
	int col = cols() + 1;
	double** matrix = AugmentedMatrix::getAugmentedMatrix(this->matrix, this->vector);
	Matrix::normalizeMatrix(matrix, row, col);
	double a8 = Matrix::maximumAbsoluteRowSumNorm(matrix, row, col - 1);
	double a1 = Matrix::maximumAbsoluteColSumNorm(matrix, row, col - 1);
	double a2 = Matrix::spectralNorm(matrix, row, col - 1);
	double q = min(3, a8, a1, a2);
	if (q > 1)
		throw QString("Cannot find solution this way. (A q > 1)");
	double* result = new double[row];
	double* remember = new double[row];
	for (int i = 0; i < row; i++) {
		result[i] = vector.get(i);
	}
	double d8 = Vector::maximumAbsoluteRowSumNorm(result, row);
	int N = (log(eps * (1 - q) / d8) / log(q)) + 1;
	for (int k = 0; k < N; k++) {
		for (int i = 0; i < row; i++) {
			//remember[i] = result[i];
			result[i] = 0;
			for (int j = 0; j < col - 1; j++)
				result[i] += matrix[i][j] * result[j];
			result[i] += matrix[i][col - 1];
			//remember[i] = result[i];
		}
	}
	Matrix::freeMatrix(matrix, row);
	Vector::freeVector(remember);
	return Vector(result, row);
}

int AugmentedMatrix::rows() const {
    return matrix.rows();
}

int AugmentedMatrix::cols() const {
    return matrix.cols();
}

double** AugmentedMatrix::getAugmentedMatrix(const Matrix& matrix, const Vector& vector) {
	int row = matrix.rows();
	int col = matrix.cols();
	double** result = new double*[row];
	for (int i = 0; i < row; i++) {
		result[i] = new double[col + 1];
		for (int j = 0; j < col; j++)
            result[i][j] = matrix(i, j);
		result[i][col] = vector.get(i);
	}
	return result;
}

double AugmentedMatrix::converge(double* v1, double* v2, int len) {
	double* array = new double[len];
	for (int i = 0; i < len; i++) {
		array[i] = (v1[i] - v2[i]);
	}
	double result = Vector::maximumAbsoluteRowSumNorm(array, len);
	Vector::freeVector(array);
	return result;
}

QString AugmentedMatrix::toString() const{
	QString result("");
	for (int i = 0; i < rows(); i++) {
		for (int j = 0; j < cols(); j++) {
            result.append(QString::number(matrix(i, j)));
			result.append(" ");
		}
		result.append("| ");
		result.append(QString::number(vector.get(i)));
		result.append("\n");
	}
	return result;
}



