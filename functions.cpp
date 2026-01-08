
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
    
    // 初始化 aux 為 rhs
    for (int l = 0; l < N; ++l) {
        aux[l] = rhs[l];
    }
    
    // 減去 A*sol
    for (int j = 1; j <= NY-1; ++j) {
        for (int i = 1; i <= NX-1; ++i) {
            const int l = ij2l(i, j, NX);
            double sum = 0.0;
            
            // 檢查每個鄰居是否存在並計算貢獻
            // 下鄰居 (i, j-1)
            if (j > 1) {
                const int l_down = ij2l(i, j-1, NX);
                sum += a10 * sol[l_down];
            }
            // 左鄰居 (i-1, j)
            if (i > 1) {
                const int l_left = ij2l(i-1, j, NX);
                sum += a01 * sol[l_left];
            }
            // 自身
            sum += a11 * sol[l];
            // 右鄰居 (i+1, j)
            if (i < NX-1) {
                const int l_right = ij2l(i+1, j, NX);
                sum += a21 * sol[l_right];
            }
            // 上鄰居 (i, j+1)
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
// 1. 預計算常數
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
    
    // 2. 初始化
    res = 2.0 * TOL;
    iters = 0;
    
    // 3. Jacobi 主迭代迴圈
    while (res > TOL && iters < MAX_ITERS) {
        // 3.1 更新所有內部點（使用 aux 暫存新解）
        for (int j = 1; j <= NY-1; ++j) {
            for (int i = 1; i <= NX-1; ++i) {
                const int l = ij2l(i, j, NX);
                double sum = 0.0;
                
                // 根據點的位置處理不同數量的鄰居
                // 這部分優化：可以分為內部區域和邊界區域來避免if判斷
                
                // 下鄰居 (i, j-1)
                if (j > 1) {
                    const int l_down = ij2l(i, j-1, NX);
                    sum += a10 * sol[l_down];
                }
                // 左鄰居 (i-1, j)
                if (i > 1) {
                    const int l_left = ij2l(i-1, j, NX);
                    sum += a01 * sol[l_left];
                }
                // 右鄰居 (i+1, j)
                if (i < NX-1) {
                    const int l_right = ij2l(i+1, j, NX);
                    sum += a21 * sol[l_right];
                }
                // 上鄰居 (i, j+1)
                if (j < NY-1) {
                    const int l_up = ij2l(i, j+1, NX);
                    sum += a12 * sol[l_up];
                }
                
                // Jacobi 更新公式
                aux[l] = (rhs[l] - sum) / a11;
            }
        }
        
        // 3.2 複製新解到 sol（同步更新）
        for (int l = 0; l < N; ++l) {
            sol[l] = aux[l];
        }
        
        // 3.3 計算殘差
        res = computeResidualMatrixFree(L, NX, NY, sol, rhs, aux);
        
        // 3.4 更新迭代計數器
        ++iters;
        
        // 可選：每100次迭代輸出進度
        if (iters % 1000 == 0) {
            LOG_FILE << "  Iteration " << iters << ", residual = " << res << endl;
        }
    }
    
    // 4. 收斂檢查與輸出
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
