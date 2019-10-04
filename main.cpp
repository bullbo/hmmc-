/**
*   @author     : Amin Bolakhrif (aminbo@kth.se)
*   @file       : main.cpp
*   @created    : today
*   @licence    : MIT
* */
#include <iostream>
#include <math.h>
#include "Matrix.hpp"
#include "HMM.hpp"

int main() {

    int rows, columns, T;

    std::cin >> rows;
    std:: cin >> columns;
    Matrix<double> A(rows, columns);
    A.assignInput();

    std::cin >> rows;
    std:: cin >> columns;
    Matrix<double> B(rows, columns);
    B.assignInput();

    std::cin >> rows;
    std:: cin >> columns;
    Matrix<double> pi(rows, columns);
    pi.assignInput();

    std:: cin >> T;
    Matrix<int> Obs(1, T);
    Obs.assignInput();
    

    HMM model(A, B, pi, Obs);
    
    double oldLogProb = INT_MIN;
    double logProb = 0;

    int maxIter = 3;


    // Iterate
    for(int i = 0; i< maxIter; i++){

        if((i < maxIter) ){  //&& (logProb > oldLogProb)){
            
            std::cout << logProb<< " " << oldLogProb << std::endl;

            oldLogProb = logProb;

            model.calcAlphaPass();
            model.calcBetaPass();
            model.calcGamma();
            model.restimateModel();

                // Compute log[P(O|lambda)]
            logProb = 0;
            for(int t = 0; t< model.T; t++){
                logProb += log(model.C.mtrx[0][t]);
            }
            
            logProb = -logProb;

        }
        }
    std::cout << model.A.toString() << std::endl;

    return 0;
}