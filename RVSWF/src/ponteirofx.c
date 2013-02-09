#include <stdio.h>
#include <stdlib.h>
#define DIMXYZ 3

void fxopera(vetargs)
  double *vetargs;
{
  int 
    i,j,k;
  double 
    *R;
    R=(double *)calloc(DIMXYZ,sizeof(double));
    for(i=0;i<DIMXYZ;i++){
      for(j=0;j<DIMXYZ;j++){
        *(R+i) = *(vetargs+j) * 2.0; 
      }
    }
    for(printf("Saida...\n"),k=0;k<DIMXYZ;k++){
      printf ("%.2f\n",*(R+k));  
    }
}

int main()
{
  double 
    *vet;
  int    
    i;
  void
    fxopera(); 
    vet=(double *)calloc(DIMXYZ,sizeof(double));
    for (i=0;i<DIMXYZ;i++){
      *(vet+i)=3.4 + 2.0;
    } 
    fxopera(vet);
}
