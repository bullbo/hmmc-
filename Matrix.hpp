#pragma once
#include <string>
#include <iostream>

template<typename T>
class Matrix {

public:

    unsigned int rows;
    unsigned int columns;
    T **mtrx;

    Matrix();
    Matrix(unsigned int m, unsigned int n);
    ~Matrix();

    /* Assign Element Values*/
    void assignInput();

    /* Operators */
    Matrix<T>& operator=(const Matrix<T>&);
    Matrix<T>& operator+=(const Matrix<T>&);
    Matrix<T>& operator-=(const Matrix<T>&);
    Matrix<T>& operator*=(const Matrix<T>&);
    Matrix<T>& operator*=(float);



    /* Returns */
    int getRows(){return rows;};
    int getColumns(){return columns;};
    std::string toString();

};

template<typename T>
Matrix<T>::Matrix(unsigned int m, unsigned int n) {
    rows = m;
    columns = n;
    mtrx = new T*[m];

    /* Initialize the matrix with 0s. */
    for(unsigned int i=0; i<rows; i++){
        mtrx[i] = new T[columns];
        for(unsigned int j=0; j<columns; j++){
            mtrx[i][j] = 0;
        }
    }
}

template<typename T>
Matrix<T>::Matrix() {
    rows = 0;
    columns = 0;
    mtrx = nullptr;
}

template<typename T>
void Matrix<T>::assignInput(){
    for(unsigned int i=0; i<rows; i++){
        for(unsigned int j=0; j< columns; j++){
            std::cin >> mtrx[i][j];
        }
    }
}

template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m) {

    rows = m.rows;
    columns = m.columns;
    mtrx = m.mtrx;
    

   return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m) {

    for (unsigned int i=0; i<rows; i++){
        for(unsigned int j=0; j<columns; j++){
            mtrx[i][j] += m.mtrx[i][j];
        }
    }
   return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m) {

    for (unsigned int i=0; i<rows; i++){
        for(unsigned int j=0; j<columns; j++){
            mtrx[i][j] -= m.mtrx[i][j];
        }
    }
   return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(float f) {

    for (unsigned int i=0; i<rows; i++){
        for(unsigned int j=0; j<columns; j++){
            mtrx[i][j] *= f;
        }
    }
   return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& m) {

    Matrix result(rows, m.columns);
    for (unsigned int i=0; i<rows; i++){
        for(unsigned int j=0; j<m.columns; j++){
            for(unsigned int k=0; k< columns; k++){
                result.mtrx[i][j] += (mtrx[i][k] * m.mtrx[k][j]);
            }
        }
    }
   return (*this = result);
}

template<typename T>
std::string Matrix<T>::toString() {
    std::string output;

    for(unsigned int i=0; i<rows; i++){
        for(unsigned int j=0; j<columns; j++){
            output += std::to_string(mtrx[i][j])+" ";
        }

        if(i !=(rows-1))
            output += "\n";
    }
    return output; 
}


template<typename T>
Matrix<T>::~Matrix() {
}


