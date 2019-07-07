#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <omp.h>
using namespace std;

/* Declearation */

#define NUM_THREAD 4
#define MAX_TIME 300000
#define RPD 4
#define k1 0.00022
#define k2 0.01875
#define k3 0.000625
#define k4 0.0625
#define k5 0.2
#define k6 0.2
#define k7 4e-8
#define k8 0.005
#define k9 0.005
#define k10 1e-4
#define k11 1e-4
#define Ka 5000
#define Kr 3500
#define OMEGA 1
#define N 2
#define P 5
#define k7_OMEGA float(k7/OMEGA)
const short int substrate[11][8]={\
{1,0,0,0,0,0,0,0},\
{0,1,0,0,0,0,0,0},\
{0,0,1,0,0,0,0,0},\
{0,0,0,1,0,0,0,0},\
{0,0,0,0,1,0,0,0},\
{0,0,0,0,0,1,0,0},\
{0,0,0,0,0,0,1,1},\
{0,0,0,0,1,0,0,0},\
{0,0,0,0,0,1,0,0},\
{0,0,0,0,0,0,1,0},\
{0,0,0,0,0,0,0,1},\
};
const short int product[11][8]={\
{1,0,0,0,1,0,0,0},\
{0,1,0,0,1,0,0,0},\
{0,0,1,0,0,1,0,0},\
{0,0,0,1,0,1,0,0},\
{0,0,0,0,1,0,1,0},\
{0,0,0,0,0,1,0,1},\
{0,0,0,0,0,0,0,1},\
{0,0,0,0,0,0,0,0},\
{0,0,0,0,0,0,0,0},\
{0,0,0,0,0,0,0,0},\
{0,0,0,0,0,0,0,0},\
};

/* Reactions & Functions */

float prng_0_1(){
    std::random_device rd;
    uniform_real_distribution<float> u(0.0,1.0);
    return u(rd);
}

void calc_a(float *spdconst,float *reactant,float *a){
    a[0] = 0.0;
    for(int n=1;n<12;n++){
        a[n] = 1.0;
    }
    for(int n=0;n<11;n++){
        for(int k=0;k<8;k++){
            if(substrate[n][k]*reactant[k] != 0.0){
                a[n+1] *= substrate[n][k]*reactant[k];
            }
        }
    }
    for(int n=1;n<12;n++){
        a[n] *= spdconst[n-1];
        a[0] += a[n];
    }
}

int reaction_sel(float *a){
    float p2 = prng_0_1();
    float sum = 0.0;
    for(int r=2;r<=11;r++){                                                 //May have problem(Solved)
        sum += a[r-1];
        if(sum <= p2*a[0] && sum+a[r] > p2*a[0]){
            return r;  }
    }
    return 1;
}

void update(float *tot_time,float *a,float *spdconst,float *reactant){
    spdconst[1] = k2*pow(reactant[7],N)/(pow(Ka,N)+pow(reactant[7],N));
    spdconst[3] = k4*pow(reactant[7],P)/(pow(Kr,P)+pow(reactant[7],P));
    calc_a(spdconst,reactant,a);
    float p1 = prng_0_1();
    float tou = (-1.0/a[0])*log(p1);
    *tot_time = *tot_time+tou;
    int sel = reaction_sel(a);
    for(int n=0;n<8;n++){
        reactant[n] = reactant[n]-substrate[sel][n]+product[sel][n];
    }
}

/* Main */

int main(){
    float reactant[8]={1e4,1e4,1e3,1e3,0,0,0,0};                            // Method of storage: PA_CONST,PR_CONST,PA,PR,mRNA_A,mRNA_R,A,R
    float a[12];                                                           
    float spdconst[11] = { k1, 0, k3, 0, k5, k6, k7_OMEGA, k8, k9, k11 };
    float tot_time=0.0;
    float(*final_data)[8] = new float[(MAX_TIME*10)][8];
    #pragma omp parallel for shared(final_data) num_threads(NUM_THREAD)
    for(int n=0;n<(MAX_TIME*10);n++){                                       //Init value to 0
        for(int k=0;k<8;k++){
            final_data[n][k] = 0;
        }
    }

    #pragma omp parallel for firstprivate(reactant,a,spdconst,tot_time) shared(final_data) num_threads(NUM_THREAD)
    for(int m=0;m<RPD;m++){                                                 //Caculation
        for(int n=0;tot_time < MAX_TIME;n++){
            int pos_rpd = (int)(tot_time*10);
            update(&tot_time,a,spdconst,reactant);
            for(int k=0;k<8;k++){
                final_data[pos_rpd][k] += reactant[k];
            }
        }
        cout << "Finished RPD " << m+1 << endl;
    }

    #pragma omp parallel for shared(final_data) num_threads(NUM_THREAD)
    for(int n=0;n<(10*MAX_TIME);n++){
        for(int p=0;p<8;p++){
            final_data[n][p] /= RPD;
        }
    }

    cout << "Writing data." << endl;

    std::ofstream out2("raw_data.txt");
    for(int n=0;n<(10*MAX_TIME);n++){
        out2 << n/10.0 << " " << final_data[n][6] << " " << final_data[n][7] << '\n';
    }
    out2.close();
    
    delete []final_data;
    return 0;
}
