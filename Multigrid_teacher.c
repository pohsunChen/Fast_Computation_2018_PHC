#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
int Jacobi_Iteration(double **U, double **F, int N);
int GaussSeidel_Iteration(double **U, double **F, int N);
int Multigrid_Iteration(double **U, double **F, int N);
double Residual(double **R, double **U, double **F, int N);
int main()
{
	int i, j, k, p, N, Max_Steps = 100, M;
	double h, **U, **F, **R, r, r0, r1;
	time_t t1, t2;

	p = 8;
	N = 1 << p;
	M = (N-1)*(N-1);
	printf("N = %d\n",N);
	printf("Total Number of Unknowns = %d\n", M);
	h = M_PI/N;
	U = (double **) malloc( (N-1)*sizeof(double*) );
	U[0] = (double *) malloc( (N-1)*(N-1)*sizeof(double) );
	for(i=1;i<N-1;++i) U[i] = U[i-1] + (N-1);
	F = (double **) malloc( (N-1)*sizeof(double*) );
	F[0] = (double *) malloc( (N-1)*(N-1)*sizeof(double) );
	for(i=1;i<N-1;++i) F[i] = F[i-1] + (N-1);
	R = (double **) malloc( (N-1)*sizeof(double*) );
	R[0] = (double *) malloc( (N-1)*(N-1)*sizeof(double) );
	for(i=1;i<N-1;++i) R[i] = R[i-1] + (N-1);

	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			U[i][j] = 0.0;
			F[i][j] = 1.0*h*h;
		}
	}
	r = r0 = Residual(R, U, F, N);
	k = 0;
	t1 = clock();
	while(r > r0*1e-6 && k < Max_Steps)
	{
		Multigrid_Iteration(U, F, N);
		r1 = Residual(R, U, F, N);
		k ++;
		printf("%d:%e %f\n",k,r1,r1/r);
		r = r1;
	}
	t2 = clock();
	printf("Total time: %f\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	return 0;
}
int Jacobi_Iteration(double **U, double **F, int N)
{
	int i, j;
	double s, **V;

	V = (double **) malloc( (N-1)*sizeof(double*) );
	V[0] = (double *) malloc( (N-1)*(N-1)*sizeof(double) );
	for(i=1;i<N-1;++i) V[i] = V[i-1] + (N-1);

	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			V[i][j] = s/(-4.0);
		}
	}
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			U[i][j] = V[i][j];
		}
	}
	free(V[0]);
	free(V);

	return 0;
}
int GaussSeidel_Iteration(double **U, double **F, int N)
{
	int i, j;
	double s;

	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			U[i][j] = s/(-4.0);
		}
	}
	for(i=0;i<N-1;++i)
	{
		for(j=N-2;j>=0;j--)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			U[i][j] = s/(-4.0);
		}
	}
	for(i=N-2;i>=0;i--)
	{
		for(j=0;j<N-1;++j)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			U[i][j] = s/(-4.0);
		}
	}
	for(i=N-2;i>=0;i--)
	{
		for(j=N-2;j>=0;j--)
		{
			s = F[i][j];
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			U[i][j] = s/(-4.0);
		}
	}
	return 0;
}
int Multigrid_Iteration(double **U, double **F, int N)
{
	int i, j;
	double s, r, **R, **Un, **Fn;
	if(N == 2)
	{
		U[0][0] = -F[0][0]/4.0;
	}
	else
	{
		// open memory for next level
		R = (double **) malloc( (N-1)*sizeof(double*) );
		R[0] = (double *) malloc( (N-1)*(N-1)*sizeof(double) );
		for(i=1;i<N-1;++i) R[i] = R[i-1] + (N-1);
		Un = (double **) malloc( (N/2-1)*sizeof(double*) );
		Un[0] = (double *) malloc( (N/2-1)*(N/2-1)*sizeof(double) );
		for(i=1;i<N/2-1;++i) Un[i] = Un[i-1] + (N/2-1);
		Fn = (double **) malloc( (N/2-1)*sizeof(double*) );
		Fn[0] = (double *) malloc( (N/2-1)*(N/2-1)*sizeof(double) );
		for(i=1;i<N/2-1;++i) Fn[i] = Fn[i-1] + (N/2-1);

		//Jacobi_Iteration(U,F,N);
		GaussSeidel_Iteration(U,F,N);
		r = Residual(R,U,F,N);
		for(i=0;i<N/2-1;++i)
		{
			for(j=0;j<N/2-1;++j)
			{
				Un[i][j] = 0.0;
				Fn[i][j] = 4*R[2*i+1][2*j+1];
				/*
				Fn[i][j] = ((R[2*i][2*j]+2*R[2*i][2*j+1]+R[2*i][2*j+2])+
				           (2*R[2*i+1][2*j]+4*R[2*i+1][2*j+1]+2*R[2*i+1][2*j+2])+
						   (R[2*i+2][2*j]+2*R[2*i+2][2*j+1]+R[2*i+2][2*j+2]))/4.0;
				*/
			}
		}
		Multigrid_Iteration(Un,Fn,N/2);
		for(i=0;i<N-1;++i)
		{
			for(j=0;j<N-1;++j)
			{
				if(i%2==1 && j%2 == 1)
				{
					R[i][j] = Un[(i-1)/2][(j-1)/2];
				}

				if(i%2==1 && j%2 == 0)
				{
					R[i][j] = 0.0;
					if(j/2-1>=0) R[i][j] += Un[(i-1)/2][j/2-1];
					if(j/2<N/2-1) R[i][j] += Un[(i-1)/2][j/2];
					R[i][j] *= 0.5;
				}
				if(i%2==0 && j%2 == 1)
				{
					R[i][j] = 0.0;
					if(i/2-1>=0) R[i][j] += Un[i/2-1][(j-1)/2];
					if(i/2<N/2-1) R[i][j] += Un[i/2][(j-1)/2];
					R[i][j] *= 0.5;
				}
				if(i%2==0 && j%2 == 0)
				{
					R[i][j] = 0.0;
					if(j/2-1>=0&&i/2-1>=0) R[i][j] += Un[i/2-1][j/2-1];
					if(j/2<N/2-1&&i/2-1>=0) R[i][j] += Un[i/2-1][j/2];
					if(j/2-1>=0&&i/2<N/2-1) R[i][j] += Un[i/2][j/2-1];
					if(j/2<N/2-1&&i/2<N/2-1) R[i][j] += Un[i/2][j/2];
					R[i][j] *= 0.25;
				}
			}
		}
		for(i=0;i<N-1;++i)
		{
			for(j=0;j<N-1;++j)
			{
				U[i][j] += R[i][j];
			}
		}
		//GaussSeidel_Iteration(U,F,N);
		free(R[0]);
		free(R);
		free(Un[0]);
		free(Un);
		free(Fn[0]);
		free(Fn);
	}
	return 0;
}
double Residual(double **R, double **U, double **F, int N)
{
	int i, j;
	double s, r=0.0;
	for(i=0;i<N-1;++i)
	{
		for(j=0;j<N-1;++j)
		{
			s = F[i][j]-U[i][j]*(-4);
			if(j-1>=0) s -= U[i][j-1];
			if(j+1< N-1) s -= U[i][j+1];
			if(i-1>=0) s -= U[i-1][j];
			if(i+1< N-1) s -= U[i+1][j];
			if(fabs(s) > r) r = fabs(s);
			R[i][j] = s;
		}
	}
	return r;
}
