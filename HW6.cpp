//Preston Peck
//CS 365
//October 30, 2017
//HW 6

#include <iostream>
#include <cmath>
using namespace std;

double N(double x);
double d1(double K, double S, double v, double r, double t0, double T);
double d2(double K, double S, double v, double r, double t0, double T);

//double europeanCallParity(double K, double S, double r, double t0, double T, double p);
double europeanPutParity(double K, double S, double r, double t0, double T, double c);
double europeanCallDividendsParity(double K, double S, double q, double r, double t0, double T, double p);
//double europeanPutDividendsParity(double K, double S, double q, double r, double t0, double T, double c);

//double europeanCall(double K, double S, double v, double r, double t0, double T);
//double europeanPut(double K, double S, double v, double r, double t0, double T);
//double europeanCallDividends(double K, double S, double v, double q, double r, double t0, double T);
//double europeanPutDividends(double K, double S, double v, double q, double r, double t0, double T);

double strike(double S, double r, double t0, double T, double c, double p);

double europeanDeltaCall(double K, double S, double v, double r, double t0, double T);
double europeanDeltaPut(double K, double S, double v, double r, double t0, double T);
double deltaFutures(double r, double t0, double T);

int main() {
    //6.1 Put-call parity
    cout << "6.1.1" << endl;
    //strike, stock, volatility, rate, time start, time end, dividend, call, put
    double K, S, v, r, t0, T, q, p, c = 0.0;
    S = 100;
    r = 10;
    K = 101;
    T = .5;
    c = 8;
    europeanPutParity(K, S, r, t0, T, c);

    cout << "6.1.2" << endl;
    S = 100;
    q = 3;
    r = 10;
    K = 101;
    T = .75;
    p = 4;
    europeanCallDividendsParity(K, S, q, r, t0, T, p);

    cout << "6.1.3" << endl;
    S = 100;
    r = 5;
    c = 6;
    p = 7;
    T = 1;
    strike(S, r, t0, T, c, p);

    cout << "6.2 Put-call parity & option pricing bounds" << endl;
    S = 20;
    r = 8;
    c = 21;
    p = 21;
    T = 1;
    strike(S, r, t0, T, c, p);
    
    cout << "6.3 Delta" << endl;
    S = 10;
    v = 50;
    r = 6;
    K = 12;
    t0 = 0.0;
    T = 0.8;
    double eDC = europeanDeltaCall(K, S, v, r, t0, T);

    K = 11;
    double eDP = europeanDeltaPut(K, S, v, r, t0, T);

    double df = deltaFutures(r, t0, T);
    Nf()
}

//HELPER Functions
double N(double x) {
  const double root = sqrt(0.5);
  return 0.5 * (1.0 + erf(x * root));
}

double d1(double K, double S, double v, double r, double t0, double T) {
    double root = sqrt(T - t0);
    return ((log(S / K) + (r * (T - t0))) / (v * root)) 
        + (.5 * v * root);
}

/*double d2(double K, double S, double v, double r, double t0, double T) {
    double root = sqrt(T - t0);
    return d1(K, S, v, r, t0, T) - (v * root);
}
*/


double europeanCallParity(double K, double S, double r, double t0, double T, double p) {
    double rate = r * .01;
    double callP = p + ((S) - (K * exp(-1 * rate * (T - t0))));

    cout << "PRICE of European CALL option: $" << callP << endl << endl;
    return callP;
}

double europeanPutParity(double K, double S, double r, double t0, double T, double c) {
    double rate = r * .01;
    double putP = c - ((S) - (K * exp(-1 * rate * (T - t0))));

    cout << "PRICE of European PUT option: $" << putP << endl << endl;
    return putP;
}

double europeanCallDividendsParity(double K, double S, double q, double r, double t0, double T, double p) {
    double rate = r * .01;
    double dividend = q * .01;

    double callPD = p + ((S * exp(-1 * dividend * (T - t0))) - (K * exp(-1 * rate * (T - t0))));
    cout << "PRICE of European CALL option w/ DIVIDENDS: $" << callPD << endl << endl;
    return callPD;
}

/*double europeanPutDividendsParity(double K, double S, double q, double r, double t0, double T, double c) {
    double rate = r * .01;
    double dividend = q * .01;

    double putPD = c - ((S * exp(-1 * dividend * (T - t0))) - (K * exp(-1 * rate * (T - t0))));
    cout << "PRICE of European PUT option w/ DIVIDENDS: $" << putPD << endl << endl;
    return putPD;
}*/

/*double europeanCall(double K, double S, double v, double r, double t0, double T) {
    double rate = r * .01;
    double volatility = v * .01;

    double call = S * N(d1(K, S, volatility, rate, t0, T)) 
        - K * exp(-1 * rate * (T - t0)) * N(d2(K, S, volatility, rate, t0, T));
    cout << "PRICE of European CALL option: $" << call << endl << endl;
    return call;
}*/

/*double europeanPut(double K, double S, double v, double r, double t0, double T) {
    double rate = r * .01;
    double volatility = v * .01;

    double put = K * exp(-1 * rate * (T - t0)) * N(-1 * d2(K, S, volatility, rate, t0, T))
        - S * N(-1 * d1(K, S, volatility, rate, t0, T));
    cout << "PRICE of European PUT option: $" << put << endl << endl;
    return put;
}*/

/*double europeanCallDividends(double K, double S, double v, double q, double r, double t0, double T) {
    double rate = r * .01;
    double dividend = q * .01;
    double volatility = v * .01;

    double callD = S * exp(-1 * dividend * (T - t0)) * N(d1(K, S, volatility, rate, t0, T)) 
        - K * exp(-1 * rate * (T - t0)) * N(d2(K, S, volatility, rate, t0, T));
    cout << "PRICE of European CALL option w/ DIVIDENDS: $" << callD << endl << endl;
    return callD;
}*/

/*double europeanPutDividends(double K, double S, double v, double q, double r, double t0, double T) {
    double rate = r * .01;
    double dividend = q * .01;
    double volatility = v * .01;

    double putD = K * exp(-1 * rate * (T - t0)) * N(-1 * d2(K, S, volatility, rate, t0, T))
        - S * exp(-1 * dividend * (T - t0)) * N(-1 * d1(K, S, volatility, rate, t0, T));
    cout << "PRICE of European PUT option w/ DIVIDENDS: $" << putD << endl << endl;
    return putD;
}*/

double strike(double S, double r, double t0, double T, double c, double p) {
    double rate = r * .01;

    double K = (-1 * c + p + S) / exp(-1 * rate * (T - t0));
    cout << "STRIKE of European options: $" << K << endl << endl;
    return K;
    
}

double europeanDeltaCall(double K, double S, double v, double r, double t0, double T) {
    double rate = r * .01;
    double volatility = v * .01;

    double deltaCall = N(d1(K, S, volatility, rate, t0, T));
    cout << "DELTA of European CALL option: " << deltaCall << endl << endl;
    return deltaCall;
}

double europeanDeltaPut(double K, double S, double v, double r, double t0, double T) {
    double rate = r * .01;
    double volatility = v * .01;

    double deltaPut = -1 * N(-1 * d1(K, S, volatility, rate, t0, T));
    cout << "DELTA of European PUT option: " << deltaPut << endl << endl;
    return deltaPut;
}

double deltaFutures(double r, double t0, double T) {
    double rate = r * .01;

    double deltaF = exp(rate * (T - t0));
    cout << "DELTA of FUTURES contract: " << deltaF << endl << endl;
    return deltaF;
}

double Nfp(double delta, double deltaF) {
    double Nf = abs((-1 * deltaCall) / deltaF);
    cout << "FUTURES contracts for DELTA ZERO PORTFOLIO: " << Nf << endl;
    return Nf;
}

double Nfp(double delta, double deltaF) {
    double Nf = abs((-1 * deltaCall) / deltaF);
    cout << "FUTURES contracts for DELTA ZERO PORTFOLIO: " << Nf << endl;
    return Nf;
}