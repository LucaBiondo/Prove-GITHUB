//
//  AnalisiMultivariata.c
//
//
//  Created by Luca Biondo on 24/04/16.
//
//


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//variabili globali

#define N 100000

//definisco i prototipi delle funzioni del programma

void acquisizione_matrice (int righe, int colonne, double M[righe][colonne]/*, char nome_file[]*/);
void prod_matrici (int righeA, int colonneA, int righeB, int colonneB, double MA[righeA][colonneA], double MB[righeB][colonneB], double MR[righeA][colonneB]);
void matrice_covarianza (int righe, int colonne, double media[colonne], double M[righe][colonne], double TAO[colonne][colonne]);
void matrice_correlazione(int righe, int colonne, double M[righe][colonne], double MCV[colonne][colonne], double MCR[colonne][colonne]);
void Matrice_Givens (int colonne, int p, int q, double G[colonne][colonne], double M[colonne][colonne]);
void Givens (int colonne, double M[colonne][colonne], double TM[colonne][colonne]);
void Householder (int colonne, double M[colonne][colonne], double TM[colonne][colonne]);
void TridiagonalizzazioneMatrice (int colonne, double M[colonne][colonne], double TM[colonne][colonne]);
void MetodoJacobi (int colonne, double M[colonne][colonne], double MJ[colonne][colonne], int I);
void MetodoQR (int colonne, double M[colonne][colonne], double MQR[colonne][colonne], int I);
void autovalori_matrice (int colonne, double M[colonne][colonne], double autovalori[colonne]);

//inizia il main

int main ()
{
    FILE *dati, *risultati;
    int righe = 0, colonne = 1, i = 0, j = 0, r = 0;
    char buffertemp[N], nome_file[N];
    dati = fopen("dati.txt", "rb");
    double temp;
    char *buff;
    char risp;
    
    printf("programma di calcolo delle matrici di covarianza e correlazione a partire da un file di dati.\n Inoltre e' possibile, tramite il calcolo degli autovalori della matrice di covarianza, calcolare il primo piano principale e rappresentare i dati proiettati su di esso.\n (metodi fino ad ora implementati: Jacobi, QR con tridiagonalizzazione della matrice tramite givens o householder)\n creatore Biondo Luca\n\n\n");
    
    //calcolo righe e colonne della matrice del file di input e creo la matrice vergine
    
    if (dati == NULL)
    {
        printf("Non posso aprire il file\n");
        exit(1);
    }
    
    fgets(buffertemp, N, dati);
    
    while (buffertemp[i] != '\n')
    {
        if (buffertemp[i] == ' ' && buffertemp[i+1] != '\n')
        {
            colonne++;
        }
        
        i++;
        
        if (i>N)
        {
            printf("il file ha troppe colonne\n");
            exit(1);
        }
    }
    
    char rbuffer[N];
    
    do
    {
        buff=fgets(rbuffer, N, dati);
        righe++;
    } while (buff!=NULL);
    
    double M[righe][colonne];
    
    printf("righe: %d, colonne: %d \n", righe, colonne);
    
    fclose (dati);
    
    //richiamo la funzione per l'acquisizione della matrice e calcolo il vettore media e le matrici di covarianza e correlazione
    
    acquisizione_matrice (righe, colonne, M);
    
    printf("matrice acquisita\n");
    
    for(i=0; i<righe; i++)
    {
        for(j=0; j<colonne; j++)
        {
            printf("%lf ", M[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    
    double TAO[colonne][colonne];
    double media[colonne];
    double MCR[colonne][colonne];
    
    
    for (i=0; i<colonne; i++)
    {
        temp=0;
        for(j=0; j<righe; j++)
        {
            temp+=M[j][i];
        }
        temp=temp/righe;
        media[i]=temp;
    }
    
    matrice_covarianza (righe, colonne, media, M, TAO);
    matrice_correlazione (righe, colonne, M, TAO, MCR);
    
    
    //stampo le matrici calcolate ed il vettore media su un file
    
    printf("vettore media\n");
    
    for(i=0; i<colonne; i++)
    {
        printf("%lf ", media[i]);
    }
    
    printf("\n");
    printf("\n");
    
    printf("matrice covarianza\n");
    
    for(i=0; i<colonne; i++)
    {
        for(j=0; j<colonne; j++)
        {
            printf("%lf ", TAO[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    
    printf("matrice correlazione\n");
    
    for(i=0; i<colonne; i++)
    {
        for(j=0; j<colonne; j++)
        {
            printf("%lf ", MCR[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    
    //richiamo l'algoritmo del calcolo degli autovalori e cerco l'autovalore maggiore
    
    printf("vuoi cercare il primo piano principale? (Y or N)\n");
    scanf("%c", &risp);
    
    if (risp=='Y')
    {
        double autovalori[colonne];
        
        autovalori_matrice (colonne, M, autovalori);
        
        printf("autovalori\n");
        
        for(i=0; i<colonne; i++)
        {
            printf("%lf ", autovalori[i]);
        }
        
        printf("\n");
        printf("\n");
        
    }else if (risp=='N')
    {
        exit(1);
    }
    
    //creo il piano di regressione e regredisco i dati
    
    //plotto il risultato su un file
  
}

void acquisizione_matrice (int righe, int colonne, double M[righe][colonne]/*, char nome_file[]*/)
{
    int i = 0, j = 0;
    FILE *dati;
    double temp;
    
    dati = fopen("dati.txt", "rb");
    
    for(i=0; i<righe; i++)
    {
        for(j=0; j<colonne; j++)
        {
            fscanf(dati, "%lf", &M[i][j]);
        }
    }
    
    fclose (dati);
}

void prod_matrici (int righeA, int colonneA, int righeB, int colonneB, double MA[righeA][colonneA], double MB[righeB][colonneB], double MR[righeA][colonneB])
{
    int i = 0, j = 0, z = 0;
    double temp1 = 0;
    
    for (i = 0; i < righeA; i++)
    {
        for (z = 0; z < colonneB; z++)
        {
            temp1=0;
            for ( j = 0; j < colonneA; j++)
            {
                temp1+=(MA[i][j]*MB[j][z]);
            }
            MR[i][z]=temp1;
        }
    }
    
}

void matrice_covarianza (int righe, int colonne, double media[colonne], double M[righe][colonne], double TAO[colonne][colonne])
{
    int i = 0, j = 0, h = 0;
    double temp = 0;
    
    for (i = 0; i < colonne; i++)
    {
        for (j = 0; j < colonne; j++)
        {
            temp=0;
            for (h = 0; h < righe; h++)
            {
                temp+=((M[h][i]*M[h][j])-(media[i]*media[j]));
            }
            temp=temp/righe;
            TAO[i][j]=temp;
        }
    }
    

    
}

void matrice_correlazione(int righe, int colonne, double M[righe][colonne], double MCV[colonne][colonne], double MCR[colonne][colonne])
{
    int i=0, j=0;
    
    for (i=0; i<colonne; i++)
    {
        for(j=0; j<colonne; j++)
        {
            MCR[i][j]=(MCV[i][j])/(sqrt(MCV[i][i]*MCV[j][j]));
        }
    }
    
    
}

void Matrice_Givens (int colonne, int p, int q, double G[colonne][colonne], double M[colonne][colonne])
{
    int i, j;
    double alpha, beta, lambda, theta;
    
    theta=(M[q][q]-M[p][p])/(2*M[p][q]);
    
    if (theta<0)
    {
        lambda= -1/(-theta+(sqrt((theta*theta)+1)));
    } else if (theta>0)
    {
        lambda= 1/(theta+(sqrt((theta*theta)+1)));
    } else
    {
        alpha=1/(sqrt(2));
        beta=1/(sqrt(2));
    }
    
    if (theta!=0)
    {
        alpha= 1/(sqrt((lambda*lambda)+1));
        beta= lambda*alpha;
    }
    
    for (i=0; i<colonne; i++)
    {
        for (j=0; j<colonne; j++)
        {
            if (i==j & i!=p & i!=q)
            {
                G[i][j]=1;
            } else if (i==p & j==p)
            {
                G[i][j]=alpha;
            } else if (i==p & j==q)
            {
                G[i][j]=beta;
            } else if (i==q & j==p)
            {
                G[i][j]=-beta;
            } else if (i==q & j==q)
            {
                G[i][j]=alpha;
            } else
            {
                G[i][j]=0;
            }
        }
    }
}

void Givens (int colonne, double M[colonne][colonne], double TM[colonne][colonne])
{
    double G[colonne][colonne];
    int i=0, j=0;
    
    for (i=1; i<colonne; i++)
    {
        for (j=i+1; j<colonne; j++)
        {
            Matrice_Givens (colonne, i, j, G, M);
            prod_matrici (colonne, colonne, colonne, colonne, G, M, TM);
            prod_matrici (colonne, colonne, colonne, colonne, TM, G, TM);
        }
    }
}

void Householder (int colonne, double M[colonne][colonne], double TM[colonne][colonne])
{
    
}

void TridiagonalizzazioneMatrice (int colonne, double M[colonne][colonne], double TM[colonne][colonne])
{
    int scelta;
    
    printf("vuoi tridiagonalizzare la matrice tramite le trasformazioni di Givens (1) o Householder(2)?\n");
    scanf("%d", &scelta);
    
    if (scelta==1)
    {
        Givens (colonne, M, TM);
    } else if (scelta==2)
    {
        Householder (colonne, M, TM);
    }
}

void MetodoJacobi (int colonne, double M[colonne][colonne], double MJ[colonne][colonne], int I)
{
    
}

void MetodoQR (int colonne, double M[colonne][colonne], double MQR[colonne][colonne], int I)
{
    
}

void autovalori_matrice (int colonne, double M[colonne][colonne], double autovalori[colonne])
{
    int trid;
    int rispp, I, i=0;
    double MA[colonne][colonne];
    
    printf("vuoi usare il metodo di Jacobi (1) o QR (2)?\n");
    scanf("%d", &rispp);
    
    printf("vuoi tridiagonalizzare la matrice dei dati prima? (1 si o 2 no)\n");
    scanf("%d", &trid);
    
    printf("quante iterazioni del metodo vuoi fare?\n");
    scanf("%d", &I);
    
    if (trid==1)
    {
        double TM[colonne][colonne];
        
        TridiagonalizzazioneMatrice (colonne, M, TM);
        
        if (rispp==1)
        {
            MetodoJacobi (colonne, TM, MA, I);
        } else if (rispp==2)
        {
            MetodoQR (colonne, TM, MA, I);
        }
    }
    
    if (rispp==1)
    {
        MetodoJacobi (colonne, M, MA, I);
    } else if (rispp==2)
    {
        MetodoQR (colonne, M, MA, I);
    }
    
    for (i=0; i<colonne; i++)
    {
        autovalori[i]=MA[i][i];
    }
    
}











