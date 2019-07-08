#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <time.h>
using namespace std;

/* Declearation */

#define NUM_THREAD 6
#define MAX_TIME 80
#define RPD 1
#define e 2.7182818
#define TOL_RANGE 0.03
#define INIT_TAU 10
#define NUM_REACTION 2
#define NUM_REACTANT 3
#define BETA 0.001
#define GAMMA 0.1

const short int substrate[NUM_REACTION][NUM_REACTANT]={\
{1,1,0},\
{0,1,0},\
};
const short int product[NUM_REACTION][NUM_REACTANT]={\
{0,2,0},\
{0,0,1},\
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
    for(int n=1;n<=NUM_REACTION;n++){
        a[n] = 1.0;
    }
    for(int n=0;n<NUM_REACTION;n++){
        for(int k=0;k<NUM_REACTANT;k++){
            if(substrate[n][k]*reactant[k] != 0.0){
                a[n+1] *= substrate[n][k]*reactant[k];
            }
        }
    }
    for(int n=1;n<=NUM_REACTION;n++){
        a[n] *= spdconst[n-1];
        a[0] += a[n];
    }
}

void calc_lambda(float *a,float delt_t,float *lambda){
    for(int n=1;n<=NUM_REACTION;n++){
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
    calc_a(spdconst,reactant,a);
    calc_lambda(a,delt_t,lambda);
    for(int reac=0;reac<NUM_REACTION;reac++){                                               // Update all the reactions
        int reac_num = calc_num(lambda[reac+1]);                                            // Lambda_0 is unavaliable (a_0 is the sum of ai)
        for(int n=0;n<NUM_REACTANT;n++){                                                    // Reactions continue for reac_num times
            reactant[n] = reactant[n]+reac_num*(-1.0*substrate[reac][n]+product[reac][n]);
            if(reactant[n]<0.0)reactant[n]=0.0;
        }
    }
}

float sel_tau(float *spdconst,float *reactant){
    float tau = 1e7;
    float temp_a[(NUM_REACTION+1)];
    float temp_reactant[NUM_REACTANT]={500,1,0};
    float diff[(NUM_REACTION+1)][NUM_REACTANT];
    calc_a(spdconst,temp_reactant,temp_a);
    for(tau=INIT_TAU;tau>0.0;tau-=1e-10){
        for(int sel=1;sel<=NUM_REACTION;sel++){
            for(int n=0;n<NUM_REACTANT;n++){
                diff[sel][n] = abs(tau*temp_a[sel]*(-1.0*substrate[sel-1][n]+product[sel-1][n]));
            }
        }
        float max = diff[1][0]/temp_a[0];
        for(int n=1;n<=NUM_REACTION;n++){
            for(int m=0;m<NUM_REACTANT;m++){
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
    float reactant[NUM_REACTANT]={500,1,0};                                                                 // Method of storage: S,I,R
    float a[(NUM_REACTION+1)];
    float lambda[(NUM_REACTION+1)];                                                                         // Lambda_0 is UNAVAILABLE                                                       
    float spdconst[NUM_REACTION] = {BETA,GAMMA};
    float delt_t = sel_tau(spdconst,reactant);                                                              // Using initial number of reactants
    float (*final_data)[NUM_REACTANT] = new float[(int)((MAX_TIME/delt_t)+1)][NUM_REACTANT];
    clock_t start, finish;
    #pragma omp parallel for shared(final_data) num_threads(NUM_THREAD)
    for(int n=0;n<(MAX_TIME*10);n++){
        for(int k=0;k<NUM_REACTANT;k++){
            final_data[n][k] = 0.0;
        }
    }

    #pragma omp parallel for firstprivate(reactant,a,spdconst,lambda) shared(final_data,delt_t) num_threads(NUM_THREAD)
    for(int m=0;m<RPD;m++){                                                                                 // Caculation of the whole system
        start = clock();
        for(int n=0;n<(int)((MAX_TIME/delt_t)+1);n++){
            update(spdconst,reactant,a,lambda,delt_t);
            for(int k=0;k<NUM_REACTANT;k++){
                final_data[n][k] += reactant[k];
            }
        }
        finish = clock();
        cout << "Finished RPD " << m+1 << " with time " << (float)(finish-start)/(float)(CLOCKS_PER_SEC) << " sec." << endl;
    }

    #pragma omp parallel for shared(final_data) num_threads(NUM_THREAD)
    for(int n=0;n<(int)((MAX_TIME/delt_t)+1);n++){
        for(int p=0;p<NUM_REACTANT;p++){
            final_data[n][p] = final_data[n][p]/RPD;
        }
    }

    cout << "Simulating system using time interval " << delt_t << endl;
    cout << "Writing data." << endl;

    std::ofstream out("SIR_tau.txt");
    for(int n=0;n<(int)((MAX_TIME/delt_t)+1);n++){
        out << n*delt_t << " " << final_data[n][0] << " " << final_data[n][1] << " " << final_data[n][2] <<  '\n';
    }
    out.close();
    
    delete []final_data;
    return 0;
}
