#ifndef MATRIX_H
#define MATRIX_H
#include <QString>
#include <QStringList>
#include <QStringBuilder>
#include <QDebug>
#include <math.h>

class Vector;

class Matrix
{
    double determinant;
    double** matrix;
    int row;
    int col;
public:
    Matrix();
    Matrix(QString text);
    Matrix(double** matrix, int row, int col);
    Matrix(double* xs, double* ys, int length);
    Matrix(const Matrix& other);
    ~Matrix();

    int rows() const;
    int cols() const;
    int footPrint() const;
    double get(int i, int j) const;
    double getDeterminant() const;
    double getMinor(int i, int j) const;
    double** getMatrix() const;

    void set(int i, int j, double value);
    double& set(int i, int j);

    Matrix inversibleMatrix() const;
    Matrix gaussianInversibleMatrix() const;
    Vector characteristicEquotation() const;
    Vector characteristicEquotation2() const;

    const Matrix multiply(double number) const;
    const Matrix multiply(const Matrix& that) const;
    const Matrix sum(const Matrix& that) const;
    const Vector multiply(const Vector& vect) const;

    Matrix& operator =(const Matrix& that);
    const Matrix operator *(const Matrix& that) const;
    const Matrix operator *(double number) const;
    const Vector operator *(const Vector& that) const;
    const Matrix operator +(const Matrix& matrix) const;
    const double operator ()(int i, int j) const;

    QString toString() const;

    static Matrix getTranslation(double a1, double a2);
    static Matrix getRotation(double angle);
    static Matrix getRotation(double x, double y, double angle);
    static Matrix getSymmetric(double x, double y);
    static Matrix getHomogeneous(double x, double y, double k);

    static double** getIdentity(int size);

    static void swapRow(double** matrix, int row1, int row2, int col);
    static void swapCol(double** matrix, int col1, int col2, int row);
    static void swapVectCol(double** matrix, const Vector& vect, int col);

    static void normalizeMatrix(double** matrix, int row, int col);
    static double maximumAbsoluteRowSumNorm(double** matrix, int row, int col);// || A ||8
    static double maximumAbsoluteColSumNorm(double** matrix, int row, int col);// || A ||1
    static double spectralNorm(double** matrix, int row, int col);// || A ||2

    static double** arrayCopy(double** array, int row, int col);

    static double getDeterminant(double **matrix, int row, int col);
    static double forwardElimination(double** augMatrix, int row, int col);
    static void backwardElimination(double** augMatrix, int row, int col);

    static double getMinor(double **matrix, int row, int col, int i, int j);
    static double** arrayCopyWithEscape(double **matrix, int row, int col, int i, int j);

    static double maxInCol(double **matrix, int rows, int colToCheck);
    static double minInCol(double **matrix, int rows, int colToCheck);

    static void freeMatrix(double** array, int row);

    static QString arrayToString(double** matrix, int row, int col);
};

#endif // MATRIX_H
