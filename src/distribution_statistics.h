#ifndef DMULTINOM_H_
#define DMULTINOM_H_

double dmultinom(int size, unsigned int *n, double *p);

unsigned int factorial(unsigned int a);

void rdirichlet( double* alpha, int size, double *theta );

/* probability density function, Dirichlet distribution */
double ddirichlet(double* x, double* alpha, int K);

#endif /*DMULTINOM_H_*/
