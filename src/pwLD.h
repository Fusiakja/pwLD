

void estimateHaploFreq( double genoFreq[4][4], double haploFreq[3][3], double *tol, int *digits, double pooHat[3] );

void whichSolutionInterval(double *poo, int *idx, double *af, double tol);

void alleleFreq( double genoFreq[4][4], double af[2]);

void whichSolutionML( double genoFreq[4][4], double *poo, int solNumb, int *validSolNumbIdx, double *af, double *Lscore);

void estimateLD(double tab[3][3], char **what, int *numb, double *LD);

void estimateGenoFreq( char **genoA, char **genoB, char **code, int *N, char **paradigm,  double Dir[9], double genoFreq[4][4],  double genoCounts[4][4], int *Neff );

void LDall(char **data, int *nrow, int *ncol, char **LD, int *LDnumb, char **code, char **paradigm, double Dir[9], double *MAF, double *tol, int *digits, double *LDmatPtr);

void bootstrapGenoFreq( double genoFreq[4][4], int *N, double genoFreqBS[4][4] );

void confidenceInterval( double genoFreq[4][4], int *N, int *nSim, char **LD, int *LDnumb, double *tol, int *digits, double *LDdist );

void credibleInterval(double genoCounts[4][4], int *nSim, double Dir[9], char **LD, int *LDnumb, double *tol, int *digits, double *LDdist);

void sampleDirichletGenoFreq( double genoCounts[4][4], double *Dir, double genoFreqRS[4][4] );

/* ###########################
	nine different likelihood terms with 
	respect to nine genotypic frequencies 
*/
double term00(double Poo);
double term01(double Poo, double pox);
double term02(double Poo, double pox);
double term10(double Poo, double pxo);
double term11(double Poo, double pox, double pxo);
double term12(double Poo, double pox, double pxo);
double term20(double Poo, double pxo);
double term21(double Poo, double pox, double pxo);
double term22(double Poo, double pox, double pxo);
