// selscan-numeric.h -- the numerical routines norm and gene used from GSL.
//
// selscan used GSL only for these, so GSL is no longer a dependency and the
// vendored lib/ and include/gsl are not needed to build.
//
// sort, sort2, quantile_sorted and mean follow the GSL 2.8 code (heapsort
// in sort/sort_source.c and sort/sort2_source.c, statistics/quantiles_source.c,
// statistics/mean_source.c), so results are unchanged. The heapsort is
// kept, not replaced by std::sort, because sort2 is not stable: the order
// in which it leaves tied keys decides the order of the second array.
//
// fit_line is ordinary least squares of y on [1, x] in centred closed form.
// GSL used an SVD (gsl_multifit_linear); estimates agree to rounding.
//
// Rng is std::mt19937_64 (its output sequence is fixed by the C++
// standard) with an unbiased bounded draw, so a --perm-seed run gives the
// same permutations on every platform. It replaces gsl_rng_default and
// gsl_ran_shuffle, so p-values for a given seed differ from GSL builds.

#ifndef __SELSCAN_NUMERIC_H__
#define __SELSCAN_NUMERIC_H__

#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

namespace snum {

inline void downheap(double *data, const size_t N, size_t k) {
    double v = data[k];
    while (k <= N / 2) {
        size_t j = 2 * k;
        if (j < N && data[j] < data[j + 1]) j++;
        if (!(v < data[j])) break;   // also stops on NaN
        data[k] = data[j];
        k = j;
    }
    data[k] = v;
}

// Ascending heapsort, as gsl_sort(data, 1, n).
inline void sort(double *data, const size_t n) {
    if (n == 0) return;
    size_t N = n - 1;
    size_t k = N / 2;
    k++;
    do {
        k--;
        downheap(data, N, k);
    } while (k > 0);
    while (N > 0) {
        std::swap(data[0], data[N]);
        N--;
        downheap(data, N, 0);
    }
}

inline void downheap2(double *data1, double *data2, const size_t N, size_t k) {
    double v1 = data1[k];
    double v2 = data2[k];
    while (k <= N / 2) {
        size_t j = 2 * k;
        if (j < N && data1[j] < data1[j + 1]) j++;
        if (!(v1 < data1[j])) break;
        data1[k] = data1[j];
        data2[k] = data2[j];
        k = j;
    }
    data1[k] = v1;
    data2[k] = v2;
}

// Sorts data1 ascending and applies the same permutation to data2, as
// gsl_sort2(data1, 1, data2, 1, n).
inline void sort2(double *data1, double *data2, const size_t n) {
    if (n == 0) return;
    size_t N = n - 1;
    size_t k = N / 2;
    k++;
    do {
        k--;
        downheap2(data1, data2, N, k);
    } while (k > 0);
    while (N > 0) {
        std::swap(data1[0], data1[N]);
        std::swap(data2[0], data2[N]);
        N--;
        downheap2(data1, data2, N, 0);
    }
}

// Quantile f of sorted data by linear interpolation at f*(n-1), as
// gsl_stats_quantile_from_sorted_data(sorted, 1, n, f). 0 when n == 0.
inline double quantile_sorted(const double *sorted, const size_t n, const double f) {
    if (n == 0) return 0.0;
    const double index = f * (n - 1);
    const size_t lhs = (int)index;
    const double delta = index - lhs;
    if (lhs == n - 1) return sorted[lhs];
    return (1 - delta) * sorted[lhs] + delta * sorted[lhs + 1];
}

// Mean by the recurrence GSL uses (long double accumulator), as
// gsl_stats_mean(data, 1, n).
inline double mean(const double *data, const size_t n) {
    long double m = 0;
    for (size_t i = 0; i < n; i++) m += (data[i] - m) / (i + 1);
    return (double)m;
}

// Least squares y = b0 + b1 x. Returns {0, 0} when x has no spread or
// fewer than 2 points.
inline std::pair<double, double> fit_line(const std::vector<double> &x, const std::vector<double> &y) {
    const size_t n = x.size();
    if (n < 2 || y.size() != n) return {0.0, 0.0};
    double mx = 0, my = 0;
    for (size_t i = 0; i < n; i++) { mx += x[i]; my += y[i]; }
    mx /= n;
    my /= n;
    double sxx = 0, sxy = 0;
    for (size_t i = 0; i < n; i++) {
        const double dx = x[i] - mx;
        sxx += dx * dx;
        sxy += dx * (y[i] - my);
    }
    if (sxx == 0) return {my, 0.0};
    const double b1 = sxy / sxx;
    return {my - b1 * mx, b1};
}

class Rng {
public:
    explicit Rng(uint64_t seed) : g(seed) {}
    // Uniform integer in [0, n), n > 0, by rejection (no modulo bias).
    uint64_t below(uint64_t n) {
        const uint64_t lim = UINT64_MAX - UINT64_MAX % n;   // a multiple of n
        uint64_t v;
        do { v = g(); } while (v >= lim);
        return v % n;
    }
    // Fisher-Yates, in the order gsl_ran_shuffle uses.
    template <class T> void shuffle(T *a, size_t n) {
        if (n < 2) return;
        for (size_t i = n - 1; i > 0; i--) std::swap(a[i], a[below(i + 1)]);
    }
private:
    std::mt19937_64 g;
};

}  // namespace snum

#endif
