#pragma once
#include <string>
#include "Matrix.hpp"
#include <vector>


class HMM {

public:
    Matrix<double> A;  //State transition Matrix.
    Matrix<double> B;  //Observation Matrix.
    Matrix<double> pi; //Initial state distribution.
    Matrix<int> O; //Initial state distribution.
    Matrix<double> C;

    Matrix<double> alpha;
    Matrix<double> beta;
    Matrix<double> delta;
    Matrix<int> deltaIdx;
    Matrix<double> gamma;
    std::vector<Matrix<double>> digamma;

    int T; //Length of the observation sequence
    int N; //Number of states in the model
    int M; //Number of observation symbols

    HMM(Matrix<double> A, Matrix<double> B, Matrix<double> pi, Matrix<int> O);
    ~HMM();

    void calcAlphaPass();
    void calcBetaPass();
    void calcGamma();

    void restimateModel();

    Matrix<int> calcViterbi();




    double probObservationSeq();

    std::string alphaToString();
    std::string betaToString();
    std::string gammaToString();


};