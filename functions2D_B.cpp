#include "functions.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

inline int ij2l(const int i, const int j, const int Nx) {
   return (j-1)*(Nx-1)+(i-1);
}

static double computeResidualMatrixFree(const double L, const int NX, const int NY,
                                       const double* const sol, const double* const rhs,
                                       double* const aux) {
    const double DX = L / static_cast<double>(NX);
    const double DY = L / static_cast<double>(NY);
    const double DXSQ = DX * DX;
    const double DYSQ = DY * DY;
    
    const double a10 = 1.0 / DYSQ;
    const double a01 = 1.0 / DXSQ;
    const double a11 = -(2.0/DXSQ + 2.0/DYSQ);
    const double a21 = 1.0 / DXSQ;
    const double a12 = 1.0 / DYSQ;
    
    double res = 0.0;
    const int N = (NX-1)*(NY-1);
    
    // Initialize aux to rhs
    for (int l = 0; l < N; ++l) {
        aux[l] = rhs[l];
    }
    
    // Subtract A*sol
    for (int j = 1; j <= NY-1; ++j) {
        for (int i = 1; i <= NX-1; ++i) {
            const int l = ij2l(i, j, NX);
            double sum = 0.0;
            
            // Check each neighbor and compute contribution
            // Bottom neighbor (i, j-1)
            if (j > 1) {
                const int l_down = ij2l(i, j-1, NX);
                sum += a10 * sol[l_down];
            }
            // Left neighbor (i-1, j)
            if (i > 1) {
                const int l_left = ij2l(i-1, j, NX);
                sum += a01 * sol[l_left];
            }
            // Self
            sum += a11 * sol[l];
            // Right neighbor (i+1, j)
            if (i < NX-1) {
                const int l_right = ij2l(i+1, j, NX);
                sum += a21 * sol[l_right];
            }
            // Top neighbor (i, j+1)
            if (j < NY-1) {
                const int l_up = ij2l(i, j+1, NX);
                sum += a12 * sol[l_up];
            }
            
            aux[l] -= sum;
            res += fabs(aux[l]);
        }
    }
    
    return res / sqrt(static_cast<double>(N));
}

int solveJacobi2D_B(const int nnz
                   ,const int N
                   ,const int * const cooRow
                   ,      int * const csrRow
                   ,const int * const cooCol
                   ,const double * const cooMat
                   ,const double* const rhs
                   ,const double TOL,const int MAX_ITERS
                   ,double* const sol
                   ,double& res , int& iters
                   ,double* const aux
                   ,std::ofstream& LOG_FILE) {

    // Preprocess to find the starting index for each row's entries in COO
    for (int i = 0; i < N; ++i) {
        csrRow[i] = -1;
    }
    int inz = 0;
    while (inz < nnz) {
        int row = cooRow[inz];
        csrRow[row] = inz;
        while (inz < nnz && cooRow[inz] == row) {
            ++inz;
        }
    }

    res = 2.0 * TOL;
    iters = 0;

    while (res > TOL && iters < MAX_ITERS) {
        // Update aux with Jacobi step using precomputed start indices
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            double a11 = 0.0;
            int p = csrRow[i];
            while (p < nnz && cooRow[p] == i) {
                int col = cooCol[p];
                double coeff = cooMat[p];
                if (col == i) {
                    a11 = coeff;
                } else {
                    sum += coeff * sol[col];
                }
                ++p;
            }
            aux[i] = (rhs[i] - sum) / a11;
        }

        // Compute residual: ||rhs - A * aux||_1 / sqrt(N)
        double temp_res = 0.0;
        for (int l = 0; l < N; ++l) {
            sol[l] = 0.0;  // Reuse sol as temp for A * aux
        }
        for (int i = 0; i < N; ++i) {
            int p = csrRow[i];
            while (p < nnz && cooRow[p] == i) {
                int col = cooCol[p];
                sol[i] += cooMat[p] * aux[col];
                ++p;
            }
        }
        for (int l = 0; l < N; ++l) {
            temp_res += fabs(rhs[l] - sol[l]);
        }
        res = temp_res / sqrt(static_cast<double>(N));

        // Copy aux back to sol
        for (int l = 0; l < N; ++l) {
            sol[l] = aux[l];
        }

        ++iters;
        if (iters % 1000 == 0) {
            LOG_FILE << "  Iteration " << iters << ", residual = " << res << endl;
        }
    }

    if (res <= TOL) {
        LOG_FILE << "Successful convergence" << endl;
        LOG_FILE << " - residual achieved " << res << endl;
        LOG_FILE << " - iterations " << iters << endl;
        return EXIT_SUCCESS;
    } else {
        LOG_FILE << "Failure of convergence procedure" << endl;
        LOG_FILE << " - residual achieved " << res << endl;
        LOG_FILE << " - iterations " << iters << endl;
        return EXIT_FAILURE;
    }
}
