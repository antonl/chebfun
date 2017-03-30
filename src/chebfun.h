/// Interface for the sample library module.
///
/// @file
#pragma once

#include "common.h"
#include <bitset>
#include <vector>
#include <tuple>
#include <functional>
#include <cmath>
#include <limits>


NAMESPACE_BEGIN(chebfun)

/// Class to represent a 1D interpolation on the interval [-1, 1] to machine
/// precision using Chebyshev polynomials. Based on Battles and Treffan (2004) paper
template <typename T> struct basic_chebfun
{
    // interpolate to this floating point type precision
    using value_t = T;
    using domain_t = std::tuple<value_t, value_t>;

    // set of flags, not sure if necessary
    std::bitset<sizeof(uint8_t)> flags;

    // domain of the original function
    domain_t domain;

    // scale of original function
    value_t scale;

    // array of chebyshev interpolant on the interval [-1, 1] (xj)
    std::vector<value_t> points;

    // array of interpolating coefficients (fj)
    std::vector<value_t> values;

    // array of interpolation weights (wj)
    std::vector<value_t> weights;

    // create the identity function f(x)=x on the default interval
    basic_chebfun();

    // create the identity function f(x)=x on an arbitrary interval [a, b]
    basic_chebfun(domain_t domain);

    // construct chebfun at interpolation points on interval given by domain
    basic_chebfun(std::vector<value_t> values, domain_t domain);

    // construct chebfun on interval in a given domain by evaluating a callable within that domain
    basic_chebfun(std::function<value_t(value_t)>, domain_t domain);

    /// Evaluate interpolant at the current point x
    value_t operator()(value_t x);


    /// Calculate the order of polynomial used
    size_t size() const;
};

/// computes the Chebyshev polynomial coefficients, with the convention that the highest coefficients are *last*
/// This is opposite to the chebfun convention, but matches the convention for using Horner's rule
template <typename value_t> std::vector<value_t> chebpoly(std::vector<value_t>& p);

NAMESPACE_BEGIN(detail)
const double pi = std::acos(-1);

// Contains flags for chebfun
enum struct Flags: uint8_t {
    Initalized,
    Prune,
};

// maps a value in domain [a, b] -> unit interval
template <typename T> inline T domain2unit(const T v, const std::tuple<T, T> domain) {
    T a, b;
    std::tie(a,b) = domain;
    return (a + b - 2.0*v)/(a-b);
}

// maps a value in unit interval [-1, 1] to domain
template <typename T> inline T unit2domain(const T v, const std::tuple<T, T> domain) {
    T a, b;
    std::tie(a,b) = domain;
    return 0.5*(b-a)*v + 0.5*(a+b);
}


template <typename T> std::vector<T> chebyshev_points(const size_t N)
{
    std::vector<T> points(N);
    for(size_t n = 0; n < N; ++n)
        points[n] = std::cos(pi*n/(N-1));

    return points;
}

template <typename T> inline std::tuple<T, T> swapdomainp(std::tuple<T, T> domain_)
{
    if(std::get<0>(domain_) > std::get<1>(domain_))
        return std::make_tuple(std::get<1>(domain_), std::get<0>(domain_));

    return domain_;
}

/// compare two doubles, example taken from cppreference
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, int ulp = 1)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) < std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
           // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}
NAMESPACE_END(detail)

NAMESPACE_END(chebfun)

