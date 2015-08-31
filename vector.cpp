#include "vector.h"

Vector::Vector(QString text) {
	if (text.size() == 0)
		throw QString("Vector text representation is empty.");
	QStringList lines = text.split("\n");
	N = lines.size();
	vector = new double[N];
	for (int i = 0; i < N; i++) {
		vector[i] = lines.at(i).toDouble();
	}
}

Vector::Vector(const Vector& that) {
	this->N = that.N;
	this->vector = that.vector;
}

Vector::Vector(double* v, int n) : vector(v), N(n) {
}

Vector& Vector::operator=(const Vector& that) {
	if (this != &that) {
		Vector::freeVector(vector);
		this->N = that.N;
		this->vector = Vector::vectorCopy(that.vector, N);
	}
	return *this;
}

Vector::~Vector() {
	Vector::freeVector(this->vector);
}

double Vector::get(int i) const {
	return vector[i];
}

int Vector::size() const {
	return N;
}

double* Vector::get() const{
	return vectorCopy(vector, N);
}

Vector Vector::operator-() const {
	double* result = new double[N];
	for (int i = 0; i < N; i++)
		result[i] = -vector[i];
	return Vector(result, N);
}

QString Vector::toString() const {
	QString result("");
	for (int i = 0; i < N; i++) {
		result.append(QString::number(vector[i], 'f'));
		result.append("\n");
	}
	return result;
}

double* Vector::vectorCopy(double *vector, int length) {
	double* result = new double[length];
	for (int i = 0; i < length; i++)
		result[i] = vector[i];
	return result;
}

QString Vector::arrayToString(double* vector, int length) {
	QString result("");
	for (int i = 0; i < length; i++) {
		result.append(QString::number(vector[i], 'f'));
		if (i < length - 1)
			result.append(", ");
	}
	return result;
}

double Vector::maximumAbsoluteRowSumNorm(double* vector, int row) {
	double result = 0;
	for (int i = 0; i < row; i++)
		if (result < fabs(vector[i]))
			result = fabs(vector[i]);
	return result;
}

void Vector::freeVector(double* vector) {
	delete[] vector;
	//vector = 0;
}
