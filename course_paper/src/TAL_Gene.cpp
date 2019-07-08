#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <time.h>
using namespace std;

/* Declearation */

#define NUM_THREAD 6
#define MAX_TIME 3000000
#define RPD 6
#define e 2.7182818
#define TOL_RANGE 0.03
#define INIT_TAU 10
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

float frac(int n){
    float result = 0.0;
    (1==n||0==n)?result=1:result=frac(n-1)*n;
    return result;
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

void calc_lambda(float *a,float delt_t,float *lambda){
    for(int n=1;n<=11;n++){
        lambda[n] = delt_t*a[n];
    }
}

int calc_num(float lambda){
    int num = 0;
    float p_sum = 0.0;
    float r = prng_0_1();
    while(p_sum <=r){
        p_sum += pow(lambda,num)*pow(e,(-1.0*lambda))/frac(num);
        num++;
    }
    return num;
}

void update(float *spdconst,float *reactant,float *a,float *lambda,float delt_t){
    spdconst[1] = k2*pow(reactant[7],N)/(pow(Ka,N)+pow(reactant[7],N));
    spdconst[3] = k4*pow(reactant[7],P)/(pow(Kr,P)+pow(reactant[7],P));
    calc_a(spdconst,reactant,a);
    calc_lambda(a,delt_t,lambda);
    for(int reac=0;reac<11;reac++){                                             // Update all the reactions
        int reac_num = calc_num(lambda[reac+1]);                                // Lambda_0 is unavaliable (a_0 is the sum of ai)
        for(int n=0;n<8;n++){                                                   // 11 Reactions continue for reac_num times
            reactant[n] = reactant[n]+reac_num*(-1.0*substrate[reac][n]+product[reac][n]);
        }
    }
}

float sel_tau(float *spdconst,float *reactant){
    float tau = 1e7;
    float temp_a[12];
    float temp_reactant[8]={1e4,1e4,1e3,1e3,0,0,0,0};
    float diff[12][8];
    calc_a(spdconst,temp_reactant,temp_a);
    for(tau=INIT_TAU;tau>0.0;tau-=1e-10){
        for(int sel=1;sel<=11;sel++){
            for(int n=0;n<8;n++){
                diff[sel][n] = abs(tau*temp_a[sel]*(-1.0*substrate[sel-1][n]+product[sel-1][n]));
            }
        }
        float max = diff[1][0]/temp_a[0];
        for(int n=1;n<=11;n++){
            for(int m=0;m<8;m++){
                if((diff[n][m]/temp_a[0])>max)max=(diff[n][m]/temp_a[0]);
            }
        }
        if(max<TOL_RANGE)break;
        tau -= 1e-3;
    }
    return tau;
}

/* Main */

int main(){
    float reactant[8]={1e4,1e4,1e3,1e3,0,0,0,0};                                // Method of storage: PA_CONST,PR_CONST,PA,PR,mRNA_A,mRNA_R,A,R
    float a[12];
    float lambda[12];                                                           // Lambda_0 is UNAVAILABLE                                                       
    float spdconst[11] = { k1, 0, k3, 0, k5, k6, k7_OMEGA, k8, k9, k11 };
    float tot_time = 0.0;
    float delt_t = sel_tau(spdconst,reactant);                                  // Using initial number of reactants
    float(*final_data)[8] = new float[(MAX_TIME*10)][8];
    cout << "Simulating system using time interval " << delt_t << endl;
    clock_t start, finish;
    #pragma omp parallel for shared(final_data) num_threads(NUM_THREAD)
    for(int n=0;n<(MAX_TIME*10);n++){
        for(int k=0;k<8;k++){
            final_data[n][k] = 0;
        }
    }

    #pragma omp parallel for firstprivate(reactant,a,spdconst,tot_time,lambda) shared(final_data,delt_t) num_threads(NUM_THREAD)
    for(int m=0;m<RPD;m++){                                                     // Caculation of the whole system
        start = clock();
        for(int n=0;tot_time < MAX_TIME;n++){
            int pos_rpd = (int)(tot_time*10);
            update(spdconst,reactant,a,lambda,delt_t);
            for(int k=0;k<8;k++){
                final_data[pos_rpd][k] += reactant[k];
            }
            tot_time += delt_t;
        }
        finish = clock();
        cout << "Finished RPD " << m+1 << " with time " << (float)(finish-start)/(double)(CLOCKS_PER_SEC) << " sec." << endl;
    }

    #pragma omp parallel for shared(final_data) num_threads(NUM_THREAD)
    for(int n=0;n<(10*MAX_TIME);n++){
        for(int p=0;p<8;p++){
            final_data[n][p] /= RPD;
        }
    }
    cout << "Writing data." << endl;

    std::ofstream out("raw_data_tau.txt");
    for(int n=0;n<(10*MAX_TIME);n++){
        out << n/10.0 << " " << final_data[n][6] << " " << final_data[n][7] << '\n';
    }
    out.close();
    
    delete []final_data;
    return 0;
}
