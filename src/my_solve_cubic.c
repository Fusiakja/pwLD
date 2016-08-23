/* poly/solve_cubic.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* solve_cubic.c - finds the real roots of x^3 + a x^2 + b x + c = 0 */

/* #include <config.h> */
#include <math.h>
#include <Rmath.h>
#include "my_solve_poly.h"
#include <R.h>
/* "-std=c99" skips this constant! */
#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif
/*
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
*/

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

int 
gsl_poly_solve_cubic (double a, double b, double c, 
                      double *x0, double *x1, double *x2, double tol, int digits)
{
	/* epsilon for comparing CR2 and CQ3. Since these values are rather big, 
		a big epsilon is used here */
	double bigtol=0.0001;

	double aRound = fprec(a, digits);
	double bRound = fprec(b, digits);
	double cRound = fprec(c, digits);
  
  double q = (aRound * aRound - 3 * bRound);
  double r = (2 * aRound * aRound * aRound - 9 * aRound * bRound + 27 * cRound);

  double Q = q / 9;
  double R = r / 54;

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  double CR2 = 729 * r * r;
  double CQ3 = 2916 * q * q * q;

	//Rprintf("a=%f, b=%f, c=%f\n", a,b,c);

	//Rprintf("\nC2R=%f\n", CR2);
	//Rprintf("\nq=%f, r=%f, Q=%f, R=%f, Q3=%f, R2=%f, CR2=%f, CQ3=%f\n",q,r,Q,R,Q3,R2,CR2,CQ3);
  
//  if (R == 0 && Q == 0)
  if( (fabs(R) < tol) && (fabs(Q) < tol) )  
  {
      *x0 = - a / 3 ;
      *x1 = - a / 3 ;
      *x2 = - a / 3 ;
      return 3 ;
    }
//  else if (CR2 == CQ3) 
    else if( fabs(CR2 - CQ3) < bigtol )
    {	//Rprintf("\ncr2 == cr3\n");
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */

      double sqrtQ = sqrt (Q);
	
      if (R > 0)
        {	
          *x0 = -2 * sqrtQ  - a / 3;
          *x1 = sqrtQ - a / 3;
          *x2 = sqrtQ - a / 3;
        }
      else
        {
          *x0 = - sqrtQ  - a / 3;
          *x1 = - sqrtQ - a / 3;
          *x2 = 2 * sqrtQ - a / 3;
        }
      return 3 ;
    }
  else if (CR2 < CQ3) /* equivalent to R2 < Q3 */
    {	//Rprintf("\ncr2 < cr3\n");
      double sqrtQ = sqrt (Q);
      double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
      double theta = acos (R / sqrtQ3);
      double norm = -2 * sqrtQ;
      *x0 = norm * cos (theta / 3) - a / 3;
      *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
      *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;
      
      /* Sort *x0, *x1, *x2 into increasing order */

      if (*x0 > *x1)
        SWAP(*x0, *x1) ;
      
      if (*x1 > *x2)
        {
          SWAP(*x1, *x2) ;
          
          if (*x0 > *x1)
            SWAP(*x0, *x1) ;
        }
      
      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
      double B = Q / A ;
      *x0 = A + B - a / 3;
      return 1;
    }
}
