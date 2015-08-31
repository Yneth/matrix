#ifndef AUGMENTEDMATRIX_H
#define AUGMENTEDMATRIX_H
#include <QString>
#include <stdarg.h>
#include "matrix.h"
#include "vector.h"

class AugmentedMatrix
{
	const Matrix& matrix;
	const Vector& vector;
public:
	//AugmentedMatrix(QString text);
	AugmentedMatrix(const AugmentedMatrix& that);
	AugmentedMatrix(const Matrix& matrix, const Vector& vector);
	~AugmentedMatrix();
	int rows() const;
	int cols() const;
	Vector cramersRule() const;
	Vector inverseMatrixRule() const;
	Vector gaussianElimination() const;
	Vector simpleIterationElimination(double eps) const;
	Vector seidelIterationElimination(double eps) const;
	QString toString() const;
	static double** getAugmentedMatrix(const Matrix& matrix, const Vector& vector);
	static double converge(double* v1, double* v2, int len);//to delete
private:
	AugmentedMatrix &operator =(const AugmentedMatrix& that);
};


inline double min(int n, ...) {
	va_list list;
	va_start(list, n);
	double min = va_arg(list, double);
	for (int i = 1; i < n; i++) {
		double value = va_arg(list, double);
		min = value < min ? value : min;
	}
	va_end(list);
	return min;
}

#endif // AUGMENTEDMATRIX_H
