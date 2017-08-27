#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <malloc.h>
#include "AB_CN.h"
#include "Allocator.h"
#include "GC.h"
#include "GC_UV.h"
#include "Linspace.h"
#include "max.h"
#include "MinMax.h"
#include "PrintV.h"
#include "PrintM.h"
#include "PrintMppp.h"
#include "PrintVV.h"
#include "Tridiagonal.h"
#include "Upwind.h"
#include "V-Cycle.h"
using namespace std;
const double pi = 4.0 * atan(1.0);
//Matrix size
const int M = 128, N =M/2, iterations=200;
double CFL = 0.6, Re = 1.0 / 0.015, Sc= 0.1,
		a=1.0, c=1.0, b=1.0, d=1.0, convergence=10E-10;

// Saving data in a temporary matrix
void matrixEq(double **vectorNew, double **vectorOld, int x, int y) {
    int x1 = x-1;
    int y1 = y-1;
    for (int i = 0; i < x-1; i++) {
        for (int j = 0; j < y-1; j++) {
            vectorOld[i][j] = vectorNew[i][j];
        }
    }
}

// Interpolation function for U and V for the outside values
void interpolation(double **vecOld, double**vecNew, int m, int n, int marker) {
	if marker == 1 {
		for (int j = 1; j < n; j++) {
			for (int i = 1; i < m; i++) {
				vecNew[i][j] = 0.5*(vecOld[i - 1][j] + vecOld[i][j]);
			}
		}
	}
	else {
		for (int j = 1; j < n; j++) {
			for (int i = 1; i < m; i++) {
				vecNew[i][j] = 0.5*(vecOld[i][j - 1] + vexOld[i][j]);
			}
		}
	}
}


int main() {	
	// Vector Definition
	vector<double> uMax, recordTimeVec, kineticEnergyVec, probeVel, Xvec, timeVec = { 6.0 };
	vector<string> names = { "t=1.txt", "t=25.txt", "t=4.txt","t=65.txt","t=8.txt","t=95.txt" };

	// Variable Initilization
	bool Flag = false;
	double xo = 0.0, xm = 2.0, yo = 0.0, ym = 1.0, omega = pi / 2.0,
		Max = 0.0, time = 0.0, Ybar,X, temp,
		dx, dy, dx2, dy2, dt_hype, t, t_u, t_v, dt, K, D1, D2, probeVel;
	int  counter = 0, counterAdj = 0, mcount1, ncount1, mcount2, ncount2, error, IM, JM1, JM2, timeVecSize = timeVec.size(),
		levels;


	int m1 = M + 1, m2 = M + 2, n1 = N + 1, n2 = N + 2,
		L = (M - 1)*N, H = (N - 1)*(M);
	levels = static_cast<int>(log(MinMax(M, N, 1)) / log(2));

	// Pointer Declaration
	double* xn = new double[m1]; double* xc = new double[m2];
	double* ynode = new double[n1]; double* yc=new double[n2];
	double* a1 = new double[L]; double* b1 = new double[L]; double* c1 = new double[L];
	double* a2 = new double[H]; double* b2 = new double[H]; double* c2 = new double[H];
	double* a3 = new double[L]; double* b3 = new double[L]; double* c3 = new double[L];
	double* a4 = new double[H]; double* b4 = new double[H]; double* c4 = new double[H];
	double* du = new double[L]; double* dv=new double[H];
	double* InfNM = new double[iterations];
	
	// x and y vectors	
		dx = (xm - xo) / M;
		dy = (ym - yo) / N;
		dx2 = dx*dx;
		dy2 = dy*dy;
		linspace(xn, dx, m1, xo, 2);
		linspace(xc, dx, m2, xo, 1);
		linspace(ynode, dy, n1, yo, 2);
		linspace(yc, dy, n2, yo, 1);
		t = 0.0;
		PrintV("xc.txt", (double *)xc, m2);
		PrintV("yc.txt", (double *)yc, n2);
		PrintV("xn.txt", (double *)xn, m1);
		PrintV("yn.txt", (double *)ynode, n1);

		// Loops to find probe
		for (int i = 0; i < m1; i++)
		{
			if (xn[i] == 1.5)
			{
				IM = i;
			}
		}
		for (int j = 0; j < n2; j++)
		{
			 if (yc[j]<0.75 && yc[j]+dy>0.75)
			{
				JM1 = j;
				JM2 = j + 1;
			}
		}

		double*** ResM = new double**[levels]; 
		double***rhs = new double**[levels];;
		double*** EPS = new double**[levels];
		double*** EPSr = new double**[levels];
		double*** EPSC = new double**[levels];
		double*** temporal = new double**[levels];
		for (int i = 0; i < levels; i++)
		{
			temporal[i] = new double*[m2];
			ResM[i] = new double*[m2];
			rhs[i] = new double*[m2];
			EPS[i] = new double*[m2];
			EPSr[i] = new double*[m2];
			EPSC[i] = new double*[m2];
			for (int j = 0; j < m2; j++)
			{
				temporal[i][j] = new double[n2];
				ResM[i][j] = new double[n2];
				rhs[i][j] = new double[n2];
				EPS[i][j] = new double[n2];
				EPSr[i][j] = new double[n2];
				EPSC[i][j] = new double[n2];
				
			}
		}
		// Required Matrices and vectors
		for (int i = 0; i < levels; i++)
		{
			for (int j = 0; j < m2; j++)
			{
				for (int k = 0; k < n2; k++)
				{
					temporal[i][j][k] = 0.0;
					ResM[i][j][k] = 0.0;
					rhs[i][j][k] = 0.0;
					EPS[i][j][k] = 0.0;
					EPSr[i][j][k] = 0.0;
					EPSC[i][j][k] = 0.0;
				}
			}
		}


		// Double Pointer Declaration and Inizialization
		double** u = new double*[m1]; double** v = new double*[m2];
		double** up = new double*[m1]; double** vp = new double*[m2];
		double** ue = new double*[m1]; double** ve = new double*[m2];
		double*** phi = new double**[1]; double*** f = new double**[1];
		double** Hun = new double*[m1]; double** Hvn = new double*[m2];
		double** Hun1 = new double*[m1]; double** Hvn1 = new double*[m2];
		double** Qu = new double*[m1]; double** Qv = new double*[m2];
		double** U = new double*[m2]; double** V = new double*[m2];
		double** magnitude = new double*[m2]; double** vorticity = new double*[m1];
		double** Y = new double*[m2]; double** dYdx = new double*[m2];
		double** dYdy = new double*[m2];
		phi[0] = new double*[m2];
		f[0] = new double*[m2];


		double** Source = new double*[m2];
		double* au = new double[H];
		double* bu = new double[H];
		double* cu = new double[H];
		double* dl = new double[H];

		for (int i = 0; i < m2; i++) {
			Source[i] = new double[n2];
		}
		for (int i = 0; i < m2; i++) {
			phi[0][i] = new double[n2];
			f[0][i] = new double[n2];
		}
		for (int i = 0; i < m2; i++) {
			v[i] = new double[n1];
			vp[i] = new double[n1];
			ve[i] = new double[n1];
			Hvn[i] = new double[n1];
			Hvn1[i] = new double[n1];
			Qv[i] = new double[n1];

			U[i] = new double[n2];
			V[i] = new double[n2];
			magnitude[i] = new double[n2];
			Y[i] = new double[n2];
			dYdx[i] = new double[n2];
			dYdy[i] = new double[n2];
		}
		for (int i = 0; i < m1; i++) {
			u[i] = new double[n2];
			up[i] = new double[n2];
			ue[i] = new double[n2];
			Hun[i] = new double[n2];
			Hun1[i] = new double[n2];
			Qu[i] = new double[n2];
			vorticity[i] = new double[n1];
		}

		for (int i = 0; i < m2; i++) {
			for (int j = 0; j < n2; j++) {
				phi[0][i][j] = 0.0;
				f[0][i][j] = 0.0;
			}
		}
		for (int i = 0; i < m2; i++) {
			for (int j = 0; j < n2; j++) {
				if (i < m1) {
					u[i][j] = 0.0;
					up[i][j] = 0.0;
					ue[i][j] = 0.0;
				}
				if (j < n1) {
					v[i][j] = 0.0;
					vp[i][j] = 0.0;
					ve[i][j] = 0.0;
				}
				if (i < m1 && j < n1) {
					vorticity[i][j] = 0.0;
				}
				magnitude[i][j] = 0.0;
				Y[i][j] = 0.0;
				U[i][j] = 0.0;
				V[i][j] = 0.0;
				dYdx[i][j] = 0.0;
				dYdy[i][j] = 0.0;
				Source[i][j] = 0.0;
			}
		}
		double xi, yi, dec1, dec2, con=0.25*0.25;

		for (int i = 0; i < m2; i++) {
			xi = xc[i];
			for (int j = 0; j < n2; j++) {
				yi = yc[j];
				dec1 = (xi - 1.5)*(xi - 1.5) + (yi - 0.5)*(yi - 0.5);
				dec2 = (xi - 0.5)*(xi - 0.5) + (yi - 0.5)*(yi - 0.5);
				if (dec1 <= con) {
					Y[i][j] = 1.0;
				}
				else if(dec2 <= con) {
					Y[i][j] = 1.0;
				}
				else {
					Y[i][j] = 0.0;
				}
			}
		}

		// Time loop start
			while (t <= 6.0) {
				ghostuv((double**)u, (double**)v, m2, n2, t, omega, (double*)xn, 1);

				// FTCS to start Adam-Bahsford
				if (t == 0.0) {
					matrixEq((double **)u, (double **)up, m1, n2);
					matrixEq((double **)v, (double **)vp, m2,n2);
				}
				counter = 0;

				// dt Definition
				t_u = MinMax(dx, dy, 1) / ( a*maxM((double **)u,m1,n2,1) + b*maxM((double **)v,m2,n1,1));
				t_v = MinMax(dx, dy, 1) / ( d*maxM((double **)v,m2,n1,1) + c*maxM((double **)u,m1,n2, 1));
				dt_hype = MinMax(t_u, t_v, 1);
				dt = CFL*dt_hype;

				// Extraction of velocities and time adjustment
				if ((t < timeVec[counterAdj]) && (t + dt > timeVec[counterAdj]) && Flag == false) {
					dt = timeVec[counterAdj] - t;
					Flag = true;
					matrixEq((double **)up, (double **)ue, m1, n2);
					matrixEq((double **)vp, (double **)ve, m1, n2);
				}
				else {
					Flag = false;
				}
				// Q & H loop
				for (int i = 0; i < m2; i++) { 			// I=0:33
					for (int j = 0; j < n2; j++) {  	// Qu
						if (i < m1) { 					// i=0:32   
							Qu[i][j] = -5.0*cos(2.1*omega*t)*(1.0 + cos(pi*MinMax(abs(8.0 * (xn[i] - 1.5)), 1.0, 1)))*(1.0 + cos(pi*MinMax(abs(16.0 * (yc[j] - 0.5)), 1.0, 1)));
							if ((i > 0 && i < m1 - 1) && (j > 0 && j < n2 - 1)) { 		// i=1:31 Hv
								Hun[i][j] = -1.0 / (4.0 * dx)*((u[i + 1][j] + u[i][j])*(u[i + 1][j] + u[i][j]) - (u[i][j] + u[i - 1][j])*(u[i][j] + u[i - 1][j])) - 1.0 / (4.0 * dy)*((u[i][j] + u[i][j + 1])*(v[i][j] + v[i + 1][j]) - (u[i][j] + u[i][j - 1])*(v[i][j - 1] + v[i + 1][j - 1]));
								Hun1[i][j] = -1.0 / (4.0 * dx)*((up[i + 1][j] + up[i][j])*(up[i + 1][j] + up[i][j]) - (up[i][j] + up[i - 1][j])*(up[i][j] + up[i - 1][j])) - 1.0 / (4.0 * dy)*((up[i][j] + up[i][j + 1])*(vp[i][j] + vp[i + 1][j]) - (up[i][j] + up[i][j - 1])*(vp[i][j - 1] + vp[i + 1][j - 1]));
							}
						}
						if (j < n1) { 					// Qv
							Qv[i][j] = 4.0*cos(1.3*omega*t)*(1.0 + cos(pi*MinMax(abs(16.0 * (xc[i] - 1.0)), 1.0, 1)))*(1.0 + cos(pi*MinMax(abs(8.0 * (ynode[j] - 0.25)), 1.0, 1)));
							if ((i>0 && i < m2 - 1) && (j>0 && j < n1 - 1)) { 			// Hv
								Hvn[i][j] = -1.0 / (4.0 * dy)*((v[i][j + 1] + v[i][j])*(v[i][j + 1] + v[i][j]) - (v[i][j] + v[i][j - 1])*(v[i][j] + v[i][j - 1])) - 1.0 / (4.0 * dx)*((u[i][j] + u[i][j + 1])*(v[i][j] + v[i + 1][j]) - (u[i - 1][j] + u[i - 1][j + 1])*(v[i - 1][j] + v[i][j]));
								Hvn1[i][j] = -1.0 / (4.0 * dy)*((vp[i][j + 1] + vp[i][j])*(vp[i][j + 1] + vp[i][j]) - (vp[i][j] + vp[i][j - 1])*(vp[i][j] + vp[i][j - 1])) - 1.0 / (4.0 * dx)*((up[i][j] + up[i][j + 1])*(vp[i][j] + vp[i + 1][j]) - (up[i - 1][j] + up[i - 1][j + 1])*(vp[i - 1][j] + vp[i][j]));
							}
						}
					}
				}

				// a, b, c and d definition
				D1 = dt / (2.0*dx2*Re);
				D2 = dt / (2.0*dy2*Re);
				for (int i = 0; i < L; i++) {
					a1[i] = -D1;
					b1[i] = 1.0 + 2.0*D1;
					c1[i] = -D1;
					du[i] = 0.0;
					a3[i] = -D2;
					b3[i] = 1.0 + 2.0*D2;
					c3[i] = -D2;
					if (i < H) {
						a2[i] = -D1;
						b2[i] = 1.0 + 2.0*D1;
						c2[i] = -D1;
						dv[i] = 0.0;
						a4[i] = -D2;
						b4[i] = 1.0 + 2.0*D2;
						c4[i] = -D2;
					}
				}
			    mcount1 = 0, ncount1 = 0, mcount2=0, ncount2=0;
				for (int i = 0; i < L; i++) {
					if (i == mcount1) {
						a1[i] = 0.0;
						c1[i + M - 2] = 0.0;
						mcount1 = mcount1 + (M - 1);
					}
					if (i == ncount1) {
						a3[i] = 0.0;
						c3[i + N - 1] = 0.0;
						b3[i] = (1 + 3 * D2);
						b3[i + N - 1] = (1 + 3 * D2);
						ncount1 = ncount1 + (N);
					}
					if (i < H) {
						if (i == mcount2) {
							a2[i] = 0.0;
							c2[i + M - 1] = 0.0;
							b2[i] = (1 + 3 * D2);
							b2[i + M - 1] = (1 + 3 * D2);
							mcount2 = mcount2 + (M);
						}
						if (i == ncount2) {
							a4[i] = 0.0;
							c4[i + N - 2] = 0.0;
							ncount2 = ncount2 + (N-1);
						}
					}
				}
	
				matrixEq((double **)u, (double **)up,m1,n2);		// Up for the next round
				matrixEq((double **)v, (double **)vp,m1,n2);		// Vp for the next round

				// ADI Method
				AB_CN((double **)u, (double **)v, (double **)Qv, (double **)Qu, (double **)Hun, (double **)Hun1, (double **)Hvn, (double **)Hvn1, (double *)a1, (double *)b1, (double *)c1, (double *)a2, (double *)b2, (double *)c2, (double *)a3, (double *)b3, (double *)c3, (double *)a4, (double *)b4, (double *)c4, dt, M, N, t, omega, (double *)xn, (double *)xc, (double *)ynode, (double *)yc, D1, D2, (double *)du, (double *)dv);
				
				// f Definition
				for (int i = 1; i < m1; i++) {
					for (int j = 1; j < n1; j++) {
						f[0][i][j] = ((u[i][j] - u[i - 1][j]) / (dx)+(v[i][j] - v[i][j - 1]) / (dy)) / dt;
					}
				}

				// V-Cycle
				VC(phi,f, M, N, dx, dy, iterations, convergence, error, (double***)ResM, (double***)rhs, (double***)EPS, (double***)EPSC, (double***)EPSr, (double***) temporal, levels, (double*) InfNM);

				// Correction of U and V using Poisson
				for (int j = 1; j < n1; j++) {
					for (int i = 1; i < m1; i++) {
						if (i < M) {
							u[i][j] = u[i][j] - dt*(phi[0][i + 1][j] - phi[0][i][j])/dx;
						}
						if (j < N) {
							v[i][j] = v[i][j] - dt*(phi[0][i][j + 1]-phi[0][i][j])/dy;
						}
					}
				}
				ghostuv((double**)u, (double**)v, m2, n2, time, omega, (double*)xn, 1);

				// Velocity Interpolation
				K = 0.0;
				time = time + dt;
				interpolation((double **)u, (double **)U,m1,n1,1);
				interpolation((double **)v, (double **)V,m1,n1,0);

				// Upwind Method
				ghostuv((double**)U, (double**)V, m2, n2, t, omega, (double*)xc, 2);
				upwind((double **)U, (double **)V, (double **)Y, (double**) dYdx, (double**) dYdy,(double**)Source, Re, Sc, dx, dy, dt, m2, n2, t,(double*)au, (double*)bu, (double*)cu, (double*)dl);
				temp = Y[1][1];
				for (int i = 1; i < m2 - 2; i++) {
					for (int j = 1; j < n2 - 2; j++) {
						temp = temp + Y[i][j];
					}
				}
				Ybar = temp / (static_cast<double>(M)*static_cast<double>(N));
				temp = 0.0;
				for (int i = 1; i < m2 - 2; i++) {
					for (int j = 1; j < n2 - 2; j++) {
						temp =  temp + (Y[i][j]-Ybar)/Ybar*(Y[i][j] - Ybar) / Ybar;
					}
				}
				X = sqrt(1.0 / (static_cast<double>(M)*static_cast<double>(N))*temp);
				Xvec.push_back(X);
				// Kinetic Energy
			    for (int i = 1; i < m1; i++) {
					for (int j = 1; j < n1; j++) {
						K = K + 0.5*(U[i][j] * U[i][j] + V[i][j] * V[i][j])*dx*dy;
					}	
				}

				// Final modifications to continue time loop after desired time modification
				if (Flag == true) {
					for (int i = 1; i < m1; i++) {
						for (int j = 1; j < n1; j++) {
							magnitude[i][j] = sqrt(U[i][j] * U[i][j] + V[i][j] * V[i][j]);
							vorticity[i][j] = (v[i + 1][j] - v[i][j]) / dx - (u[i][j + 1] - u[i][j]) / dx;
						}
					}
					t = t - dt;
					matrixEq((double **)up, (double **)u, m1, n2);
					matrixEq((double **)vp, (double **)v, m2, n1);
					matrixEq((double **)ve, (double **)vp, m2, n1);
					matrixEq((double **)ue, (double **)up, m1, n2);
					if (counterAdj < timeVecSize - 1) {
						counterAdj++;
					}
				}

				// Probe Velocity Extraction
				probeVel = 0.5*(u[IM][JM1] + u[IM][JM2]);
				probeVel.push_back(probeVel);

				// Finding Max U
				if (t >0 && t <= 2.5) { 
					Temp = probeVel;
					if (Temp > Max) {
						Max = Temp;
					}
				}
				
				recordTimeVec.push_back(t);
				kineticEnergyVec.push_back(K);
			}	// Time loop ends

			// Clearing vectors before new round
			delete[] a1; delete[] b1; delete[] c1;
			delete[] a2; delete[] b2; delete[] c2;
			delete[] a3; delete[] b3; delete[] c3;
			delete[] a4; delete[] b4; delete[] c4;
			
			// Print of Kinetic Energy vector, U vector, time and the max U 
			uMax.push_back(Max);
			PrintVV("Uprobe.txt", probeVel);
			PrintVV("K.txt", kineticEnergyVec);
			PrintVV("time.txt", recordTimeVec);
			PrintVV("UGCI.txt", uMax);
		return 0;
}