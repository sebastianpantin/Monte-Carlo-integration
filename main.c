/* C-program with functions for generation of
N(mu,sigma^2)-distributed random variables */
#include <math.h> //for mathfunctions
#include <stdio.h>
#include <stdlib.h>
#include <time.h> // for timing

double * approximateIntegral(double n){
    double S_N, errorTerm,x_i, y;
    static double answer[2];
    S_N = 0;
    errorTerm = 0;
    for(int i=0; i<n; i++){
        x_i = drand48();//normrnd(mu,sigma);
        y = 4/(1+x_i*x_i);
        S_N = S_N + y;
        errorTerm = errorTerm + y*y;
    }
    S_N=S_N/n;
    answer[1] = S_N;
    answer[2] = errorTerm;
    return answer;


}

double * importanceSampling(double n)
{
    double S_N, sigmaSum, x_i, y, p, errorTerm, variance,u;
    static double answer[3];
    S_N = 0;
    sigmaSum = 0;
    errorTerm = 0;
    for(int i=0; i<n; i++){
        u=drand48();
        x_i = 2-sqrt(4-3*u); //Inverse CDF
        y = 4/(1+x_i*x_i);
        p = (4-2*x_i)/3;
        S_N = S_N + y/p;
        sigmaSum = sigmaSum + pow((y/p),2);
    }
    S_N = S_N/n;
    variance = sigmaSum/n - pow(S_N,2);
    errorTerm = sqrt(variance/n);
    answer[1] = S_N;
    answer[2] = errorTerm;
    answer[3] = variance;
    return answer;
}

double * controlVariates(double n)
{
    double S_N, sigmaSum, x_i, f, g, variance;
    static double answer[3];
    S_N = 0;
    sigmaSum = 0;
    variance = 0;
    for(int i=0; i<n; i++){
        x_i = drand48();//normrnd(mu,sigma);
        f = 4/(1+x_i*x_i);
        g = (4-2*x_i);
        S_N = S_N + f - g;
        sigmaSum = sigmaSum + pow((f-g+3.0),2);

    }

    S_N = 3 + S_N/n;
    variance = sigmaSum/n - pow(S_N,2);
    answer[1] = S_N;
    answer[2] = sqrt(variance/n);
    answer[3] = variance;
    return answer;
}

double * antitheticVariates(double n)
{
    double S_N, sigmaSum, x_i,f, fneg, variance, errorTerm;
    static double answer[3];
    S_N = 0;
    sigmaSum = 0;
    for(int i=0; i<n; i++){
        x_i = drand48();//normrnd(mu,sigma);
        f = 4/(1+x_i*x_i);
        fneg = 4/(1+(1-x_i)*(1-x_i));
        S_N = S_N + f/2 + fneg/2;
        sigmaSum = sigmaSum + pow((f/2+fneg/2),2);
    }
    S_N = S_N/n;
    variance = sigmaSum/n - pow(S_N,2);
    errorTerm = sqrt(variance/n);
    answer[1] = S_N;
    answer[2] = errorTerm;
    answer[3] = variance ;
    return answer;
}

double * stratifiedSampling(double n)
{
    double S_N, sigmaSum, x_i, f, errorTerm, nInterval, variance;
    static double answer[3];
    S_N = 0;
    errorTerm = 0;
    nInterval = 4;
    variance = 0;
    for(int j = 1; j<= nInterval; j++) {
        double a = (j-1)/nInterval;
        double b = j/nInterval;
        double S_nIntervall = 0;
        double sigmaSumSquare = 0;
        sigmaSum = 0;
        for (int i = 0; i < n/nInterval; i++) {
            x_i = (b-a)*drand48() + a;//normrnd(mu,sigma);
            f = 4 / (1 + x_i * x_i);
            S_nIntervall = S_nIntervall + f;
            sigmaSum = sigmaSum + f;
            sigmaSumSquare = sigmaSumSquare + pow(f,2);
        }
        S_N = S_N+(b-a)/n*S_nIntervall;
        variance = variance + sigmaSumSquare/(n)-pow(sigmaSum/(n),2);
        errorTerm = errorTerm + pow((b-a),2)/n * variance;
    }
    answer[1] = S_N * 4;
    answer[2] = sqrt(errorTerm);
    answer[3] = variance;
    return answer;
}

double uniformRandom()
{
    return ( (double)(rand()) + 1.)/( (double)(RAND_MAX) + 1. );
}


//function for generate normal random variable (mu,sigma) (Box-Muller)
double normrnd(double mu,double sigma)
{
    double chi,eta,z;
    chi=drand48();
    eta=drand48();
    z=mu+sigma*sqrt(-2*log(chi))*cos(2*M_PI*eta);
    return z;
}

main()
{
// Define parameters and conduct Monte Carlo simulations
    srand48(time(NULL)); //Make a new seed for the random number generator
    double *answer;
    double piError, error;
    double n = pow(10,7);

    answer = approximateIntegral(n);
    piError  = *(answer+1)-M_PI;
    error = sqrt(*(answer+2)/n-*(answer+1)**(answer+1))/sqrt(n);
    printf("\n Standard: \n Integral: %f \n Error est: %f \n Real error: %f \n \n", *(answer + 1), error, piError);

    answer = importanceSampling(n);
    piError = M_PI - *(answer+1);
    printf("\n Importance sampling: \n Integral: %f \n Error est: %f \n Real error: %f \n \n", *(answer + 1),*(answer + 2), piError);

    answer = controlVariates(n);
    piError = M_PI - *(answer+1);
    printf("\n Control Variates: \n Integral: %f \n Error est: %f \n Real error: %f\n \n", *(answer + 1), *(answer + 2), piError);

    answer = antitheticVariates(n);
    piError = M_PI - *(answer+1);
    printf("\n Antithetic Variates: \n Integral: %f \n Error est: %f \n Real error: %f\n\n", *(answer + 1), *(answer + 2), piError);

    answer = stratifiedSampling(n);
    piError = M_PI - *(answer+1);
    printf("\n Stratified Sampling: \n Integral: %f \n Est error: %f \n Real error: %f\n\n", *(answer + 1), *(answer + 2), piError);

}
