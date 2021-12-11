#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void runCalc(float***,int,int,int,int,float,float,float);
void outputData(FILE*,float***,int,int,int);
int determineSections(int,int,int);



int main(int argc, char const *argv[])
{

    int i,j,k,T = atoi(argv[1]);

    FILE* FPTR;
    char filename[20] = "outputData.txt";

    // MxNxO matrix
    int M = atoi(argv[2]);  //i
    int N = M;              //j
    int O = atoi(argv[3]);  //k


    // utt = c^2 (uxx + uyy)
    float initial_value = 10.0; //initial displacement
    float dt = 0.1;
    float dh = 0.16;
    float c = 1.0;

    // Define 3D Array
    printf("Allocating memory\n");
    float*** A = (float***)malloc(M*sizeof(float**));

    // Allocate Memory for 3D Array
    for(i=0; i<M; i++)
    {
        A[i] = (float**)malloc(N*sizeof(float*));

        for(j=0; j<N; j++)
        {
            A[i][j] = (float*)malloc(O*sizeof(float));
        }
    }

    // Assign Initial Values to Array
    for (i=0; i<M; i++)
    {
        for (j=0; j<N; j++)
        {
            for (k=0; k<O; k++)
            {
                A[i][j][k]= 0.0;
            }
        }
    }
    // Set initial displacement at center of grid
    A[(M/2)][(N/2)][0] = initial_value;

    FPTR = fopen(filename,"w");

    char str;


    //distribute sections before entering parallel section
    int sectionWidth = (int)round(M/(T));
    int Begin[T], Stop[T];
    for (int t=0; t<T; t++)
    {
        Begin[t] = (t)*sectionWidth+1;
        Stop[t] = (t+1)*sectionWidth+1;
    }
    Stop[T-1] = M-1;  //ensure last row is correct


    // Do Calculations and Output Data
    int startRow,endRow;    //individual start row and end row for each thread number

    //output first set of data before entering parallel
    outputData(FPTR,A,M,N,0);

    //Begin timer
    double execTime, endTime;
    double startTime = omp_get_wtime();

    #pragma omp parallel num_threads(T) shared(Begin,Stop,A,M,N,dt,dh,c,T) private(startRow,endRow,k)
    {
        //assign start row and end row
        int tid = omp_get_thread_num();
        startRow = Begin[tid];
        endRow = Stop[tid];

        for (k=1; k<O; k++){
            runCalc(A,startRow,endRow,N,k,dt,dh,c);
            #pragma omp barrier
        }
    }

    //measure calculation time
    printf("Completed Calculation\n");
    endTime = omp_get_wtime();
    execTime = (endTime-startTime);

    //output data
    
    for (k=0;k<O;k++){
        outputData(FPTR,A,M,N,k);
    }
    

    fclose(FPTR);

    //print thread number, grid size, # of time steps, and real time elapsed
    printf("~~~ DETAILS ~~~\n");
    printf("Thread #:\t%d\n",T);
    printf("Grid size:\t%dx%d\n",M,N);
    printf("Exec. Time:\t %3.7f sec\n",execTime);

    // Free Memory
    for(i=0; i<M; i++){
        for(j=0; j<N; j++){
            free(A[i][j]);
        }
        free(A[i]);
    }
    free(A);


    return 0;

}

void runCalc(float*** A, int startRow, int endRow, int N, int k, float dt, float dh, float c)
{
    int i,j;

    for(i = startRow; i<(endRow); i++)
    {
        for(j=1; j<N-1; j++)
        {
            A[i][j][k+1] = 2.0*A[i][j][k]-A[i][j][k-1] + (pow(dt,2)*pow(c,2)/pow(dh,2))*(A[i+1][j][k]+A[i-1][j][k]+A[i][j+1][k]+A[i][j-1][k]-4.0*A[i][j][k]);
        }
    }
    return;
}

void outputData(FILE* FPTR, float*** A, int M, int N, int k)
{
    int i,j;
    for(i=0; i<M; i++)
    {
        for(j=0; j<N; j++)
        {
            fprintf(FPTR,"%.6f,",A[i][j][k]);
        }
        fprintf(FPTR,"\n");
    }
    fprintf(FPTR,"\n");
}







