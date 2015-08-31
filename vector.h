#ifndef VECTOR_H
#define VECTOR_H
#include <QString>
#include <QStringList>
#include <QDebug>
#include <math.h>

class Vector {
public:
    Vector() { }
	Vector(QString text);
	Vector(double* v, int n);
	Vector(const Vector &other);
	Vector& operator=(const Vector& that);
	Vector operator-() const;
	~Vector();
	double* get() const;
	double get(int i) const;
	int size() const;
	QString toString() const;
	static double maximumAbsoluteRowSumNorm(double* vector, int row); // ||A||8
	static void freeVector(double* vector);
	static double* vectorCopy(double* vector, int length);
	static QString arrayToString(double* vector, int length);
private:
	double* vector;
	int N;
};

#endif // VECTOR_H
