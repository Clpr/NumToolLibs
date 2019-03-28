#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main()
{
    int tic,toc;
    float *tmp;
    float *simu_BM( int T, int N, float mu, float sigma );


    tic = time();
    float tmp2[200 * 100];
    tmp2 = simu_BM( 200, 100, 0.0, 0.1 );
    toc = time();

    printf( "Eclipse: %d", &toc - &tic );

    return 0;
}



float *simu_BM( int T, int N, float mu, float sigma )
{
    // declare
    float DeltaT = 1 / N;
    static float Res[ T*N ];
    float pThreshold = 0.5 * ( 1.0 + mu / sigma * sqrt(DeltaT) );
    float hRise = sigma * sqrt(DeltaT);
    // initialize Res (X0 = 0)
    Res[0] = 0.0;
    // generate Delta X
    for( int x = 0; x < (T*N-1); x++ )
    {
        if( rand() <= pThreshold )
        {
            Res[x+1] = Res[x] + hRise;
        }else{
            Res[x+1] = Res[x] - hRise;
        }
    }

    return Res;
}





























