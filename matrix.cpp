#include "matrix.h"
#include "augmentedmatrix.h"

Matrix::Matrix() : matrix(0), row(0), col(0)
{

}

Matrix::Matrix(QString text) : determinant(0) {
    if (text.size() == 0)
        throw QString("Text representation is empty.");

    QRegExp lineEscape("\n");
    QStringList lines = text.split(lineEscape);
    row = lines.size();

    QRegExp spaceEscape("[^0-9.-]");
    QStringList line = lines.at(0).split(spaceEscape, QString::SkipEmptyParts);
    col = line.size();

    matrix = new double*[row];
    for (int i = 0; i < row; i++) {
        matrix[i] = new double[col];
        for (int j = 0; j < col; j++) {
            matrix[i][j] = line.at(j).toDouble();
        }
        if (i != row - 1) line = lines.at(i + 1).split(spaceEscape, QString::SkipEmptyParts);
    }
}

Matrix::Matrix(double** matrix, int row, int col) : determinant(0) {
    if (row <= 0 || col <= 0)
        throw QString("WrongBoundsException.");
    this->row = row;
    this->col = col;

    this->matrix = matrix;
}

Matrix::Matrix(const Matrix &other) {
    this->determinant = other.determinant;
    Matrix::freeMatrix(this->matrix, this->row);
    this->row = other.row;
    this->col = other.col;
    this->matrix = Matrix::arrayCopy(other.matrix, this->row, this->col);
}

Matrix::Matrix(double *xs, double *ys, int length)
{
    this->row = length;
    this->col = 2;
    this->matrix = new double*[length];
    for (int i = 0; i < length; i++) {
        this->matrix[i] = new double[2];
    }
    for (int i = 0; i < length; i++) {
        this->matrix[i][0] = xs[i];
        this->matrix[i][1] = ys[i];
    }
}

Matrix::~Matrix() {
    Matrix::freeMatrix(this->matrix, this->row);
}

Matrix& Matrix::operator =(const Matrix& that) {
    if (this == &that)
        return *this;
    Matrix::freeMatrix(this->matrix, this->row);
    this->row = that.row;
    this->col = that.col;
    this->matrix = Matrix::arrayCopy(that.matrix, row, col);
    return *this;
}

const Matrix Matrix::operator *(const Matrix& that) const {
    return this->multiply(that);
}

const Matrix Matrix::multiply(const Matrix& that) const {
    if (col != that.rows())
        throw QString("Cannot multiply this matrixes.");

    double** result = new double*[row];
    for (int i = 0; i < row; i++) {
        result[i] = new double[this->cols()];
        for (int j = 0; j < that.cols(); j++) {
            result[i][j] = 0;
            for (int k = 0; k < that.cols(); k++) {
                result[i][j] += this->get(i, k) * that.get(k, j);
            }
        }
    }
    return Matrix(result, this->rows(), that.cols());
}

const Matrix Matrix::operator *(double number) const {
    return this->multiply(number);
}

const Matrix Matrix::multiply(double number) const {
    double** result = new double*[row];
    for (int i = 0; i < row; i++) {
        result[i] = new double[col];
        for (int j = 0; j < col; j++)
            result[i][j] = matrix[i][j] * number;
    }
    return Matrix(result, row, col);
}

const Matrix Matrix::operator +(const Matrix& that) const {
    return this->sum(that);
}

const Matrix Matrix::sum(const Matrix& that) const {
    if (this->row != that.rows() || this->col != that.cols())
        throw QString("Cannot sum this matrixes");

    double** result = new double*[row];
    for (int i = 0; i < row; i++) {
        result[i] = new double[col];
        for (int j = 0; j < col; j++) {
            result[i][j] = this->get(i, j) + that.get(i, j);
        }
    }
    return Matrix(result, row, col);
}

const Vector Matrix::operator *(const Vector& vect) const {
    return this->multiply(vect);
}

const Vector Matrix::multiply(const Vector& vect) const {
    if (col < vect.size())
        throw QString("Cannot multiply");
    double* result = new double[row];
    for (int i = 0; i < row; i++) {
        result[i] = 0;
        for (int j = 0; j < col; j++)
            result[i] += matrix[i][j] * vect.get(j);
    }
    return Vector(result, row);
}

double Matrix::get(int i, int j) const {
    if (i >= row || j >= col)
        throw QString("Out Of Bounds Exception");
    return matrix[i][j];
}

const double Matrix::operator ()(int i, int j) const {
    return this->get(i, j);
}

Vector Matrix::characteristicEquotation() const {
    if (row != col)
        throw new QString("Cannot find equotation");
    double* randomVector = new double[row];
    for (int i = 0; i < row; i++)
        randomVector[i] = i == 0 ? 1 : 0;

    QVector<Vector*> vectors;
    vectors.append(new Vector(randomVector, row));

    qDebug() << vectors[0]->toString();
    double** matrix = new double*[row];
    for (int i = 0; i < row; i++) {
        matrix[i] = new double[col];
        vectors.append(new Vector(this->multiply(*vectors.at(i))));
        qDebug() << vectors[i+1]->toString();
    }

    for (int i = 0; i < row; i++){
        Vector* temp = vectors.at(vectors.size() - 2 - i);
        for (int j = 0; j < col; j++) {
            matrix[j][i] = temp->get(j);
        }
    }
    Matrix m(matrix, row, col);
    Vector right = -*vectors.at(row);
    AugmentedMatrix am(m, right);
    qDebug() << am.toString();
    return am.cramersRule();
}

Vector Matrix::characteristicEquotation2() const {
    if (row != col)
        throw new QString("Cannot find equotation");
    double* prints = new double[row];
    prints[0] = this->footPrint();
    Matrix temp(*this);
    for (int i = 1; i < row; i++) {
        temp = temp.multiply(*this); // temp = temp * (*this);
        prints[i] = temp.footPrint();
    }
    double* result = new double[row];
    for (int i = 0; i < row; i++) {
        double sum = prints[i];
        for (int j = i - 1, index = 0; j >= 0; j--, index++) {
            sum += prints[j] * result[index];
        }
        result[i] = -(1.0 / (i + 1))*sum;
    }
    return Vector(result, row);
}

int Matrix::footPrint() const {
    int print = 0;
    for (int i = 0; i < row; i++)
        print += matrix[i][i];
    return print;
}

double** Matrix::getMatrix() const {
    return Matrix::arrayCopy(matrix, row, col);
}

void Matrix::set(int i, int j, double value)
{
    if (i > rows() || j > cols() || i < 0 || j < 0) {
        throw QString("Out of bound exception.");
    }
    matrix[i][j] = value;
}

double &Matrix::set(int i, int j)
{
    return matrix[i][j];
}

double Matrix::getDeterminant() const {
    if (this->determinant != 0)
        return determinant;
    return getDeterminant(matrix, row, col);
}

double Matrix::getDeterminant(double** array, int row, int col) {
    double** matrix = Matrix::arrayCopy(array, row, col);
    double determinant = forwardElimination(matrix, row, col);
    freeMatrix(matrix, row);
    return determinant;
}

double Matrix::forwardElimination(double** matrix, int row, int col) {
    qDebug() << "---------------------";
    double determinant = 1;
    for (int i = 0; i < row; i++) {
        int big = i;
        for (int l = i; l < row; l++) {
            if (fabs(matrix[big][i]) < fabs(matrix[l][i])) {
                big = l;
            }
        }
        if (matrix[big][i] != matrix[i][i]) {
            swapRow(matrix, i, big, col);
            determinant *= -1;
        }
        double ii = matrix[i][i];
        determinant *= (ii * pow(-1, i + i));
        if (ii == 0) {
            break;
        }
        for (int j = i; j < col; j++) {
            matrix[i][j] /= ii;
        }
        for (int k = i + 1; k < row; k++) {
            double first = matrix[k][i];
            for (int j = i; j < col; j++) {
                matrix[k][j] -= (matrix[i][j] * first);
            }
        }
        qDebug() << arrayToString(matrix, row, col*2);
    }
    return determinant;
}

void Matrix::backwardElimination(double** augMatrix, int row, int col) {
    for (int i = row - 1; i > 0; i--) {
        for (int k = i - 1; k >= 0; k--) {
            double first = augMatrix[k][i];
            for (int j = i; j < col; j++) {
                augMatrix[k][j] -= (augMatrix[i][j] * first);
            }
        }
    }
}

double Matrix::getMinor(int i, int j) const {
    return getMinor(matrix, row, col, i, j);
}

double Matrix::getMinor(double **matrix, int row, int col, int i, int j) {
    double** factor = arrayCopyWithEscape(matrix, row, col, i, j);
    double determinant = getDeterminant(factor, row - 1, col - 1);
    freeMatrix(factor, row - 1);
    return determinant;
}

Matrix Matrix::inversibleMatrix() const {
    double determinant = getDeterminant();
    if (determinant == 0)
        throw QString("cannot get inversion.");
    double** factorMatrix = new double*[row];
    for (int i = 0; i < row; i++) {
        factorMatrix[i] = new double[col];
        for (int j = 0; j < col; j++)
            factorMatrix[i][j] = (pow(-1, i + j) * getMinor(j, i)) / determinant;
    }
    return Matrix(factorMatrix, row, col);
}

Matrix Matrix::gaussianInversibleMatrix() const {
    double determinant = getDeterminant();
    if (determinant == 0)
        throw QString("cannot get inversion.");
    int row = rows();
    int col = 2 * cols();
    double** augMatrix = new double*[row];
    for (int i = 0; i < row; i++) {
        augMatrix[i] = new double[col];
        for (int j = 0; j < col / 2; j++)
            augMatrix[i][j] = get(i, j);
        for (int j = col / 2; j < col; j++)
            if (i == j % cols())// cols()
                augMatrix[i][j] = 1;
            else
                augMatrix[i][j] = 0;
    }
    forwardElimination(augMatrix, row, col);
    backwardElimination(augMatrix, row, col);
    double** result = new double*[row];
    for (int i = 0; i < row; i++) {
        result[i] = new double[col / 2];
        for (int j = col / 2; j < col; j++) {
            result[i][j % cols()] = augMatrix[i][j]; // cols()
        }
    }
    freeMatrix(augMatrix, row);
    return Matrix(result, row, col / 2);
}

double** Matrix::arrayCopyWithEscape(double **matrix, int row, int col, int i, int j) {
    double** result = new double*[row - 1];
    int m = -1;
    int n = -1;
    for (int k = 0; k < row; k++) {
        if (k == i) continue;
        m++;
        result[m] = new double[col - 1];
        for (int l = 0; l < col; l++) {
            if (l == j) continue;
            n++;
            result[m][n] = matrix[k][l];
        }
        n = -1;
    }
    return result;
}

double Matrix::maxInCol(double **matrix, int rows, int colToCheck)
{
    double max = matrix[0][colToCheck];
    for (int i = 0; i < rows; i++)
        if (max < matrix[i][colToCheck])
            max = matrix[i][colToCheck];
    return max;
}

double Matrix::minInCol(double **matrix, int rows, int colToCheck)
{
    double min = matrix[0][colToCheck];
    for (int i = 0; i < rows; i++)
        if (min > matrix[i][colToCheck])
            min = matrix[i][colToCheck];
    return min;
}

QString Matrix::arrayToString(double** matrix, int row, int col) {
    QString result;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            result.append(QString::number(matrix[i][j]));
            result.append(" ");
        }
        result.append("\n");
    }
    return result;
}

void Matrix::normalizeMatrix(double** matrix, int row, int col) {
    for (int i = 0; i < row; i++) {
        double ii = matrix[i][i];
        for (int j = 0; j < col; j++)
            matrix[i][j] /= (ii * (-1));
        matrix[i][i] = 0;
    }
}

void Matrix::swapCol(double **matrix, int col1, int col2, int row) {
    for (int i = 0; i < row; i++) {
        double temp = matrix[i][col1];
        matrix[i][col1] = matrix[i][col2];
        matrix[i][col2] = temp;
    }
}

void Matrix::swapRow(double **matrix, int row1, int row2, int col) {
    for (int j = 0; j < col; j++) {
        double temp = matrix[row1][j];
        matrix[row1][j] = matrix[row2][j];
        matrix[row2][j] = temp;
    }
}

void Matrix::swapVectCol(double **matrix, const Vector& vect, int col) {
    for (int i = 0; i < vect.size(); i++)
        matrix[i][col] = vect.get(i);
}

double** Matrix::arrayCopy(double **array, int row, int col) {
    double** result = new double*[row];
    for (int i = 0; i < row; i++) {
        result[i] = new double[col];
        for (int j = 0; j < col; j++) {
            result[i][j] = array[i][j];
        }
    }
    return result;
}

double Matrix::maximumAbsoluteRowSumNorm(double** matrix, int row, int col) { // || A ||8
    double result = 0;
    for (int i = 0; i < row; i++) {
        double sum = 0;
        for (int j = 0; j < col; j++)
            sum += fabs(matrix[i][j]);
        if (sum > result)
            result = sum;
    }
    return result;
}

double Matrix::maximumAbsoluteColSumNorm(double** matrix, int row, int col) { // || A ||1
    double result = 0;
    for (int j = 0; j < col; j++) {
        double sum = 0;
        for (int i = 0; i < row; i++)
            sum += fabs(matrix[i][j]);
        if (sum > result)
            result = sum;
    }
    return result;
}

double Matrix::spectralNorm(double** matrix, int row, int col) { // || A ||2
    double result = 0;
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            result += matrix[i][j] * matrix[i][j];
    return sqrt(result);
}

int Matrix::rows() const{
    return row;
}

int Matrix::cols() const{
    return col;
}

QString Matrix::toString() const {
    QString result("");
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            result.append(QString::number(matrix[i][j]));
            result.append(" ");
        }
        result.append("\n");
    }
    return result;
}

Matrix Matrix::getTranslation(double a1, double a2)
{
    int size = 3;
    double** matrix = getIdentity(size);

    matrix[2][0] = a1;
    matrix[2][1] = a2;

    return Matrix(matrix, size, size);
}

Matrix Matrix::getRotation(double angle)
{
    int size = 3;
    double** rotation = getIdentity(size);
    double radian = angle * 3.1415926 / 180.0;
    rotation[0][0] = cos(radian);
    rotation[0][1] = sin(radian);
    rotation[1][0] = -sin(radian);
    rotation[1][1] = cos(radian);

    return Matrix(rotation, size, size);
}

Matrix Matrix::getRotation(double x, double y, double angle)
{
    Matrix rotation(getTranslation(-x, -y));
    rotation = rotation.multiply(getRotation(angle));
    rotation = rotation.multiply(getTranslation(x, y));

    return rotation;
}

Matrix Matrix::getSymmetric(double x, double y)
{
    int size = 3;
    double** matrix = getIdentity(size);

    matrix[0][0] = -1;
    matrix[1][1] = -1;
    matrix[2][0] = 2*x;
    matrix[2][1] = 2*y;

    return Matrix(matrix, size, size);
}

Matrix Matrix::getHomogeneous(double x, double y, double k)
{
    int size = 3;
    double** matrix = getIdentity(size);

    matrix[0][0] = k;
    matrix[1][1] = k;
    matrix[2][0] = (1.0-k)*x;
    matrix[2][0] = (1.0-k)*y;

    return Matrix(matrix, size, size);
}

double** Matrix::getIdentity(int size)
{
    double** matrix = new double*[size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new double[size];
        std::fill_n(matrix[i], size, 0);
        matrix[i][i] = 1;
    }
    return matrix;
}

void Matrix::freeMatrix(double** array, int row) {
    if (array == 0 || row <= 1) {
        return;
    }
    for (int i = 0; i < row; i++) {
        if (array[i] == 0) {
            continue;
        }
        delete[] array[i];
        array[i] = 0;
    }
    delete[] array;
    array = 0;
}
