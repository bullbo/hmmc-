#include "HMM.hpp"
#include "Matrix.hpp"
#include <iostream>
using namespace std;

HMM::HMM(Matrix<double> A, Matrix<double> B, Matrix<double> pi, Matrix<int> O) {

    this->A  = A;
    this->B  = B;
    this->pi  = pi;
    this->O  = O;
    
    T = O.getColumns();
    N = pi.getColumns();
    M = B.getColumns();

    alpha = Matrix<double>(T, N);
    beta = Matrix<double>(T, N);
    gamma = Matrix<double>(T,N);
    delta = Matrix<double>(T, N);
    deltaIdx = Matrix<int>(T, N);

    C = Matrix<double>(1, T);

    // Initialize digamma vector
    for(int t=0; t<T-1; t++){
        Matrix<double> temp(N, N);
        digamma.push_back(temp);
    }
}

void HMM::calcAlphaPass(){
    // Compute alpha_0(i), the probability of the partial observation sequence up
    // to time 0, where the underlying Makrov process is in state i at t=0.
    for(int i=0; i < N; i++){
        C.mtrx[0][0] = 0;
        int Obs_0 = O.mtrx[0][0];
        
        alpha.mtrx[0][i] = pi.mtrx[0][i]*B.mtrx[i][Obs_0];
        C.mtrx[0][0] += alpha.mtrx[0][i];
    }

    // Scaling alpha_0(i) to prevent tendency to 0 as T increases.  
    C.mtrx[0][0] = 1/C.mtrx[0][0];
    for(int i= 0; i< N; i++){
        alpha.mtrx[0][i] *= C.mtrx[0][0];
    }
    
    // Compute alpha_t(i)
    for(int t=1; t<T; t++){
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                alpha.mtrx[t][i] += alpha.mtrx[t-1][j] * A.mtrx[j][i];
            }
            int Obs_t = O.mtrx[0][t];
            alpha.mtrx[t][i] *= B.mtrx[i][Obs_t];
            C.mtrx[0][t]+= alpha.mtrx[t][i];
        }

        //Scaling alpha_t(i).
        C.mtrx[0][t] = 1/C.mtrx[0][t];
        for(int i=0; i<N; i++){
            alpha.mtrx[t][i] *=C.mtrx[0][t];
        }
        
    }
}


double HMM::probObservationSeq(){
    double result = 0.0;
    for(int i=0; i<N; i++){
        result += alpha.mtrx[T-1][i];
    }
    return result;
}

void HMM::calcBetaPass(){

    // Let Beta_{T-1}(i) = 1, scaled by c_{T-1}
    for(int i=0; i< N-1; i++){
        beta.mtrx[T-1][i] = C.mtrx[0][T-1];
    }

    // Beta-pass
    for(int t= T-2; t>= 0; t--){
        for(int i=0; i < N-1; i++){
            for(int j = 0; j < N-1 ; j++){
                beta.mtrx[t][i] += A.mtrx[i][j]*B.mtrx[j][t+1]*beta.mtrx[t+1][j];
            }
            // Scale Beta elements with same scale factor as the corresponding Alpha element
            beta.mtrx[t][i] = C.mtrx[0][t]*beta.mtrx[t][i];
        }
    }

}

Matrix<int> HMM::calcViterbi(){
    for(int i=0; i<N; i++){
        delta.mtrx[0][i] = B.mtrx[i][O.mtrx[0][0]]*pi.mtrx[0][i];
    }

    double currMax = 0, prevMax = 0;

    for(int t= 1; t<T; t++){
        for(int i=0; i<N; i++){
            prevMax = 0;
            for(int j = 0; j< N; j++){
                currMax = A.mtrx[j][i]*delta.mtrx[t-1][j]*B.mtrx[i][O.mtrx[0][t]];
                if(currMax >= prevMax){
                    prevMax = currMax;
                    delta.mtrx[t][i] = currMax;
                    deltaIdx.mtrx[t][i] = j;

                }
            }
        }
    }

    double probMax = 0, rowMax = 0;

    for(int i=0; i<N; i++){
        if(delta.mtrx[T-1][i] > probMax){
            probMax = delta.mtrx[T-1][i];
            rowMax = i;
        }
    }

    Matrix<int> seq(1, T);
    int tempRowMax = rowMax;

    for(int t=T-1; t>= 0; t--){
        seq.mtrx[0][t] = tempRowMax;
        tempRowMax = deltaIdx.mtrx[t][tempRowMax];
    }

    return seq;
}

void HMM::calcGamma(){

    double denom;
    for(int t=0; t<T-1; t++){
        denom = 0;
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                denom += alpha.mtrx[t][i]*A.mtrx[i][j]*B.mtrx[j][O.mtrx[0][t+1]]*beta.mtrx[t+1][j];
            }
        }

        for(int i=0; i<N; i++){
            for(int j=0; j< N; j++){
                digamma[t].mtrx[i][j] = alpha.mtrx[t][i]*A.mtrx[i][j]*B.mtrx[j][O.mtrx[0][t+1]]*beta.mtrx[t+1][j];
                digamma[t].mtrx[i][j] /= denom;
                gamma.mtrx[t][i] += digamma[t].mtrx[i][j];
            }
        }
    }

    // Special case for gamma_{T-1}(i)
    denom = 0;
    for(int i = 0; i<N; i++){
        denom += alpha.mtrx[T-1][i];
    }

    if(denom == 0){denom = 1;}

    for(int i=0; i<N; i++){
        gamma.mtrx[T-1][i] = alpha.mtrx[T-1][i]/denom;
    }
}

void HMM::restimateModel(){

    double numer, denom = 0;

    // Re-estimate pi
    for(int i=0; i<N; i++){
        pi.mtrx[0][i] = gamma.mtrx[0][i];
    }

    // Re-estimate A
    for(int i=0; i< N; i++){
        for(int j=0; j<N; j++){
            numer = 0;
            denom = 0;
            for(int t=0; t < T-1; t++){
                numer += digamma[t].mtrx[i][j];
                denom += gamma.mtrx[t][i];
            }
            A.mtrx[i][j] = numer/denom;
        }
    }

    // Re-estimate B
    for(int i = 0; i< N; i++){
        for(int j=0; j<M; j++){
            numer = 0;
            denom = 0;
            for(int t = 0; t< T; t++){
                if(O.mtrx[0][t] == j){
                    numer += gamma.mtrx[t][i];
                }
                denom += gamma.mtrx[t][i];
            }
            B.mtrx[i][j] = numer/denom;
        }
    }
}


std::string HMM::betaToString(){
    return beta.toString();
}

std::string HMM::gammaToString(){
    return gamma.toString();
}

HMM::~HMM() {
}