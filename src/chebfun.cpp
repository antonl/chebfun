/// Implementation of the sample library module.
///
#include "chebfun.h"

NAMESPACE_BEGIN(chebfun)

template <typename T> basic_chebfun<T>::basic_chebfun()
    : flags(0), domain{-1.0, 1.0}, scale(1.0), points{-1.0, 1.0}, values{-1.0, 1.0}, weights{0.5, -0.5}
{
}

template <typename value_t> basic_chebfun<value_t>::basic_chebfun(basic_chebfun::domain_t domain_)
    : flags(0), points{-1.0, 1.0}
{
    // check that domain is in correct order
    domain = detail::swapdomainp(domain_);
    value_t a, b;
    std::tie(a, b) = domain;
    weights = {0.5, -0.5};
    values = {a, b};
    scale = std::max({std::abs(a), std::abs(b)});
}

template <typename value_t> basic_chebfun<value_t>::basic_chebfun(std::vector<value_t> values_, basic_chebfun::domain_t domain_)
    : flags(0)
{
    // check that domain is in correct order
    domain = detail::swapdomainp(domain_);

    std::for_each(std::begin(values_), std::end(values_),
                  [this](value_t v){scale = std::abs(v) > scale? std::abs(v) : scale;});

    values.resize(values_.size());
    std::transform(std::begin(values_), std::end(values_), std::begin(values), [this](value_t v){return v/scale;});

    points = detail::chebyshev_points<value_t>(values.size());

    weights.resize(values.size());
    std::fill(std::begin(weights), std::end(weights), 1.0);
    weights.front() = weights.back() = 0.5;

    // advance iterator by 2, flipping the sign as we go
    for (auto it = std::next(std::begin(weights)); it < std::end(weights); std::next(it, 2))
        (*it) *= -1;
}

template <typename T> basic_chebfun<T>::basic_chebfun(std::function<value_t (value_t)>, basic_chebfun::domain_t domain_)
{
    // check that domain is in correct order
    domain = detail::swapdomainp(domain_);
}

/// implement barycentric interpolation formula
template <typename value_t> value_t basic_chebfun<value_t>::operator()(value_t x) {
    // transform point to unit interval
    auto ux = detail::domain2unit(x, domain);

    // check that the value we're trying to calculate is not one of the interpolations
    //auto it = std::find_if(std::begin(points), std::end(points),
    //                       [&ux](value_t v){return detail::almost_equal(ux, v);});
    auto it = std::find(std::begin(points), std::end(points), ux); // apparently a less complicated approach also works

    if (it != std::end(points)) { // exact match
        auto index = std::distance(std::begin(points), it);
        return values[index]; // we're done
    }

    // would be nice to have zip functionality to iterate over weights, points, and values at the same time
    // have to use loop instead ;(
    // this will crash if all vectors are not the same size, btw

    // variables to accumulate the denominator
    value_t num, denom, diff;
    num = denom = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        diff = ux - points[i];
        diff = weights[i]/diff;
        denom += diff;
        num += diff*values[i];
    }
    return num/denom;
}

template basic_chebfun<double>::basic_chebfun();
template basic_chebfun<double>::basic_chebfun(std::tuple<double, double>);
template basic_chebfun<double>::basic_chebfun(std::vector<double>, std::tuple<double, double>);
template double basic_chebfun<double>::operator()(double);
template basic_chebfun<float>::basic_chebfun();
NAMESPACE_END(chebfun)
