#include <iostream>
#include <random>
#include "functions.h"

using namespace std;

inline int ij2l(const int i, const int j, const int Nx) {
   return (j-1)*(Nx-1)+(i-1);
}

int solveJacobi2D_C(const double L
                   ,const int NX, const int NY
                   ,const double TOL,const int MAX_ITERS
                   ,double* const sol
                   ,const double* const rhs
                   ,double& res , int& iters
                   ,double* const aux
                   ,std::ofstream& LOG_FILE) {

   for (int i = 1; i < NX; ++i) {          

         
         sol[ij2l(i,1,NX)] = 0.0;                // south boundary
         sol[ij2l(1,i,NX)] = 0.0;                // west boundary
         sol[ij2l(NX-1,i,NX)] = 0.0;             // east boundary
         sol[ij2l(i,NX-1,NX)] = sin(M_PI*i*L/NX);// north boundary

      }

   res = 2*TOL;
   const int N = (NX-1)*(NY-1);

   iters = 0;
   while (res>TOL) {
         
      
      for (int j=2; j<NY-1; ++j) {
         for (int i=2; i<NX-1; ++i) {
         const int l11 = ij2l(i,j,NX);
         const int l01 = ij2l(i,j-1,NX);
         const int l10 = ij2l(i-1,j,NX);
         const int l21 = ij2l(i+1,j,NX);
         const int l12 = ij2l(i,j+1,NX);

         aux[l11] =(sol[l21] + sol[l10] + sol[l12] + sol[l01]) * 0.25;
         }
      }


      res = 0.0;            
      // Compute the diffrence between the current and last solution. Not the residual
/*      
      for (int l = 0; l<N;++l) {
      res += pow(sol[l]-aux[l],2);
      }
      res= sqrt(res/N);
*/

      for (int j=2; j<NY-1; ++j) {
         for (int i=2; i<NX-1; ++i) {
         const int l11 = ij2l(i,j,NX);
         const int l01 = ij2l(i,j-1,NX);
         const int l10 = ij2l(i-1,j,NX);
         const int l21 = ij2l(i+1,j,NX);
         const int l12 = ij2l(i,j+1,NX);

         res += abs(aux[l21] + aux[l10] + aux[l12] + aux[l01] - 4*aux[l11]);
         }
      }
      res /= sqrt(N);

      for (int l = 0; l<N;++l) {
         sol[l] = aux[l];
      }

      ++iters;
      if (iters == MAX_ITERS) break;

   }
   if (iters < MAX_ITERS) {
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

