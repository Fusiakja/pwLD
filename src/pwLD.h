void estimateFrequencies(  double genoFreq[3][3],  double haploFreq[3][3], int *Neff,  double Dir[4], int *mc, char **paradigm);

void estimateHaploFreq(  double genoFreq[4][4],  double haploFreq[3][3],  double *tol, int *digits,  double pooHat[3] );

void whichSolutionInterval( double *poo, int *idx,  double *af,  double tol);

void alleleFreq(  double genoFreq[4][4],  double af[2]);

void whichSolutionML(  double genoFreq[4][4],  double *poo, int solNumb, int *validSolNumbIdx,  double *af,  double *Lscore);

void estimateLD( double tab[3][3], char **what, int *numb,   double *LD,  double *HSweight, double *exponent);

void estimateGenoFreq( char **genoA, char **genoB, char **code, int *N,  char **paradigm,  double Dir[4],  double genoFreq[4][4],  double genoCounts[4][4], int *Neff);

void LDall(  char **data, int *nrow, int *ncol, char **LD, int *LDnumb,  char **code, char **paradigm,  double *Dir,  double *MAF,  double *tol, int *digits,  double *LDmatPtr,  double *HSweight, int *ci, int *mc, char **strategy,  double *alpha,  double *cilow,  double *ciup, int *nsim, double *LDdist, double *vars, int *intervall);

void bootstrapHapFreq(  double genoFreq[3][3], int *N,  double genoFreqBS[3][3], char **paradigm );

void bootstrapGenoFreq(  double genoFreq[4][4], int *N,  double genoFreqBS[4][4], char **paradigm );

void confidenceInterval(  double genoCounts[3][3], double genoFreq[3][3], int *N, char **paradigm,int *nSim, char **LD, int *LDnumb,  double *tol, int *digits,   double *LDdist,  double *HSweight,  double *alpha, char **strategy,   double *cilow,   double *ciup,  double tables[*nSim][4],  double Dir[4],   double vars[*LDnumb],int *intervall, int *mc);

void confidenceGenoInterval(  double genoCounts[4][4], double genoFreq[4][4], int *N, char **paradigm,int *nSim, char **LD, int *LDnumb,  double *tol, int *digits,  double *LDdisti,  double *HSweight,  double *alpha, char **strategy,  double *cilow,  double *ciup,  double tables[*nSim][9],  double Dir[4],  double vars[*LDnumb],int *intervall, int *mc);

void credibleInterval(int genoCounts[3][3], int *nSim,  double Dir[4], char **LD, int *LDnumb,  double *tol, int *digits,  double *LDdist,  double *HSweight,  double *alpha,  double *cilow,  double *ciup);

void sampleDirichletGenoFreq( int genoCounts[3][3],  double *Dir,  double genoFreqRS[3][3] );

double quantile(int sims, double LDs[], double alpha);

void sort(  double numbers[], int count);

void MIG(  char **data, int *nrow, int *ncol, char **LD, int *LDnumb, char **code, char **paradigm,  double *Dir,  double *MAF,  double *tol, int *digits,  double *LDmatPtr,  double *HSweight,  double *ci, int *mc, char **strategy,  double *alpha,  double *cilow,  double *ciup, int *nsim,  double *LDdist,  double *MIG1,  double *MIG2);

/* ###########################
 nine different likelihood terms with 
 respect to nine genotypic frequencies 
 */
double term00( double Poo);
double term01( double Poo,  double pox);
double term02( double Poo,  double pox);
double term10( double Poo,  double pxo);
double term11( double Poo,  double pox,  double pxo);
double term12( double Poo,  double pox,  double pxo);
double term20( double Poo,  double pxo);
double term21( double Poo,  double pox,  double pxo);
double term22( double Poo,  double pox,  double pxo);
