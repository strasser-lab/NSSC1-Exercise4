
#include "functions.h"
#include <cmath>     // New add
#include <iostream>  // New add

using namespace std;

inline int ij2l(const int i, const int j, const int Nx) {
   return (j-1)*(Nx-1)+(i-1);
}

//
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
    
    // Initialize aux as rhs
    for (int l = 0; l < N; ++l) {
        aux[l] = rhs[l];
    }
    
    // Subtract A*sol
    for (int j = 1; j <= NY-1; ++j) {
        for (int i = 1; i <= NX-1; ++i) {
            const int l = ij2l(i, j, NX);
            double sum = 0.0;
            
            // Check each neighbor points is existence and compute contribution
            // Bottom neighbor point (i, j-1)
            if (j > 1) {
                const int l_down = ij2l(i, j-1, NX);
                sum += a10 * sol[l_down];
            }
            // Left neighbor point (i-1, j)
            if (i > 1) {
                const int l_left = ij2l(i-1, j, NX);
                sum += a01 * sol[l_left];
            }
            // Self point
            sum += a11 * sol[l];
            // Right neighbor point (i+1, j)
            if (i < NX-1) {
                const int l_right = ij2l(i+1, j, NX);
                sum += a21 * sol[l_right];
            }
            // Top neighbor point (i, j+1)
            if (j < NY-1) {
                const int l_up = ij2l(i, j+1, NX);
                sum += a12 * sol[l_up];
            }
            
            aux[l] -= sum;
            res += aux[l] * aux[l];
        }
    }
    
    return sqrt(res / N);
}
//

int solveJacobi2D_C(const double L
                   ,const int NX, const int NY
                   ,const double TOL,const int MAX_ITERS
                   ,double* const sol
                   ,const double* const rhs
                   ,double& res , int& iters
                   ,double* const aux
                   ,std::ofstream& LOG_FILE) {
//
// 1. Precompute constants
    const double DX = L / static_cast<double>(NX);
    const double DY = L / static_cast<double>(NY);
    const double DXSQ = DX * DX;
    const double DYSQ = DY * DY;
    
    const double a10 = 1.0 / DYSQ;
    const double a01 = 1.0 / DXSQ;
    const double a11 = -(2.0/DXSQ + 2.0/DYSQ);
    const double a21 = 1.0 / DXSQ;
    const double a12 = 1.0 / DYSQ;
    
    const int N = (NX-1)*(NY-1);
    
    // 2. Initialization
    res = 2.0 * TOL;
    iters = 0;
    
    // 3. Jacobi main iteration loop
    while (res > TOL && iters < MAX_ITERS) {
        // 3.1 Update all interior points (use aux to store new solution temporarily)
        for (int j = 1; j <= NY-1; ++j) {
            for (int i = 1; i <= NX-1; ++i) {
                const int l = ij2l(i, j, NX);
                double sum = 0.0;
                
                // Handle different number of neighbor points based on point location
                // Optimization：can separate interior and boundary regions to avoid if statements
                
                // Bottom neighbor point (i, j-1)
                if (j > 1) {
                    const int l_down = ij2l(i, j-1, NX);
                    sum += a10 * sol[l_down];
                }
                // Left neighbor point (i-1, j)
                if (i > 1) {
                    const int l_left = ij2l(i-1, j, NX);
                    sum += a01 * sol[l_left];
                }
                // Right neighbor point (i+1, j)
                if (i < NX-1) {
                    const int l_right = ij2l(i+1, j, NX);
                    sum += a21 * sol[l_right];
                }
                // Top neighbor point (i, j+1)
                if (j < NY-1) {
                    const int l_up = ij2l(i, j+1, NX);
                    sum += a12 * sol[l_up];
                }
                
                // Jacobi formula update
                aux[l] = (rhs[l] - sum) / a11;
            }
        }
        
        // 3.2 Copy new solution to sol (synchronous update)
        for (int l = 0; l < N; ++l) {
            sol[l] = aux[l];
        }
        
        // 3.3 Compute residual
        res = computeResidualMatrixFree(L, NX, NY, sol, rhs, aux);
        
        // 3.4 Update iteration counter
        ++iters;
        
        // output progress every 1000 iterations
        if (iters % 1000 == 0) {
            LOG_FILE << "  Iteration " << iters << ", residual = " << res << endl;
        }
    }
    
    // 4. Convergence check and output
    if (res <= TOL) {


      ofstream SOL("Sol.log");

      for (int j=1; j<NY; ++j) {
          for (int i=1; i<NX; ++i) {
            if (i != NX-1){
              SOL << sol[ij2l(i,j,NX)] << ",";
            }
            else{
              SOL << sol[ij2l(i,j,NX)]; 
            }
          }     

         if (j == NY-1){

         }
         else{
           SOL << endl;
         }
       
      }

        LOG_FILE << "Successfull convergence" << endl;
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
