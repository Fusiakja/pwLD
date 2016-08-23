/*
 *  Copyright (C) 2006  Hin-Tak Leung
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available via WWW at
 *  http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
 *  writing to the Free Software Foundation, Inc., 51 Franklin Street
 *  Fifth Floor, Boston, MA 02110-1301  USA.
 */

/* 
   The two routines were taken from GSL 1.7, with minor changes
   for isolated operation without the rest of GNU Scientific Library. 
   This header file is a replacement header to hook up the GSL code 
   with the rest.
*/

int
gsl_poly_solve_cubic (double a, double b, double c,
                      double *x0, double *x1, double *x2, double tol, int digits);

int
gsl_poly_solve_quadratic (double a, double b, double c,
                          double *x0, double *x1);
