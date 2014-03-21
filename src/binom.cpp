#include "binom.h"

long double fact(int x)
{
    return floor(0.5 + exp(factln(x)));
}

//From Numerical Recipes in C
//returns n choose k
long double nCk(int n, int k)
{
    return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)));
}

//From Numerical Recipes in C
//returns ln(gamma(xx)) for xx >0
long double gammln(long double xx)
{
    long double x, y, tmp, ser;
    static long double cof[6] = {76.18009172947146, -86.50532032941677,
                                 24.01409824083091, -1.231739572450155,
                                 0.1208650973866179e-2, -0.5395239384953e-5
                                };
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

//From Numerical Recipes in C
//returns ln(n!)
long double factln(int n)
{
    static long double a[101];

    if (n <= 1) return 0.0;
    if (n <= 100) return a[n] ? a[n] : (a[n] = gammln(n + 1.0));
    else return gammln(n + 1.0);
}

