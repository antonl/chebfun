/// Test suite for the sample module.
///
/// Link all test files with the `gtest_main` library to create a command-line 
/// test runner.
///
#include <gtest/gtest.h>
#include "chebfun.h"

using namespace chebfun;

TEST(Utils, MapsUnit2DomainAndBack) {
    auto domainf = std::make_tuple(0.f, 1.f);
    ASSERT_FLOAT_EQ(0.5f, detail::unit2domain(0.0f, domainf));
    ASSERT_FLOAT_EQ(0.0f, detail::unit2domain(-1.0f, domainf));
    ASSERT_FLOAT_EQ(1.0f, detail::unit2domain(1.0f, domainf));

    ASSERT_FLOAT_EQ(-1.0f, detail::domain2unit(0.0f, domainf));
    ASSERT_FLOAT_EQ(0.0f, detail::domain2unit(0.5f, domainf));
    ASSERT_FLOAT_EQ(1.0f, detail::domain2unit(1.0f, domainf));

    auto domain = std::make_tuple(0.0, 1.0);
    ASSERT_DOUBLE_EQ(0.5, detail::unit2domain(0.0, domain));
    ASSERT_DOUBLE_EQ(0.0, detail::unit2domain(-1.0, domain));
    ASSERT_DOUBLE_EQ(1.0, detail::unit2domain(1.0, domain));

    ASSERT_DOUBLE_EQ(-1.0, detail::domain2unit(0.0, domain));
    ASSERT_DOUBLE_EQ(0.0, detail::domain2unit(0.5, domain));
    ASSERT_DOUBLE_EQ(1.0, detail::domain2unit(1.0, domain));
}

TEST(ChebfunBasic, DoubleConstructorInvariants) {
    double a, b;

    auto cb1 = basic_chebfun<double>();
    ASSERT_EQ(cb1.values.size(), 2);
    ASSERT_EQ(cb1.values.size(), cb1.points.size());
    ASSERT_DOUBLE_EQ(cb1.scale, 1.0);

    auto cb2 = basic_chebfun<double>({10.0, -5.0});
    ASSERT_EQ(cb2.values.size(), 2);
    ASSERT_EQ(cb2.values.size(), cb2.points.size());
    ASSERT_DOUBLE_EQ(cb2.scale, 10.0);
    std::tie(a, b) = cb2.domain;
    ASSERT_GT(b, a);

    auto cb3 = basic_chebfun<double>({-5.0, -10.});
    std::tie(a, b) = cb3.domain;
    ASSERT_GT(b, a);
    ASSERT_DOUBLE_EQ(-10.0, a);

    auto cb4 = basic_chebfun<double> ({1}, {-5, 10});
    ASSERT_DOUBLE_EQ(cb4.scale, 1.0);
    std::tie(a, b) = cb4.domain;
    ASSERT_DOUBLE_EQ(cb4.scale, 1.0);
    ASSERT_EQ(cb4.values.size(), 1);
    ASSERT_EQ(cb4.values.size(), cb4.points.size());

}

TEST(ChebfunBasic, IdentityInterpolation) {
    auto abstol = 10*std::numeric_limits<double>::epsilon();
    auto cb1 = basic_chebfun<double>();

    ASSERT_DOUBLE_EQ(cb1(0.), 0.);
    ASSERT_DOUBLE_EQ(cb1(0.5), 0.5);
    ASSERT_DOUBLE_EQ(cb1(-0.5), -0.5);

    // check fast path when x is one of the interpolation points
    ASSERT_DOUBLE_EQ(cb1(1.0), 1.0);
    ASSERT_DOUBLE_EQ(cb1(-1.0), -1.0);

    auto cb2 = basic_chebfun<double>({-5., 10.});
    //ASSERT_TRUE(detail::almost_equal(cb2(0.0), 0., 4));
    ASSERT_TRUE(std::abs(cb2(0.1)- 0.1) < abstol);
    ASSERT_TRUE(std::abs(cb2(0.5)- 0.5) < abstol);
    ASSERT_TRUE(std::abs(cb2(-0.5)+ 0.5) < abstol);

    // check boundaries
    ASSERT_TRUE(std::abs(cb2(-5.)+ 5) < abstol);
    ASSERT_TRUE(std::abs(cb2(10)- 10) < abstol);
}

TEST(ChebfunUtils, ChebyshevPoints) {
    auto pts1 = detail::chebyshev_points<double>(7);
    ASSERT_DOUBLE_EQ(pts1.front(), -pts1.back());
    ASSERT_LT(std::abs(pts1[3]), std::numeric_limits<double>::epsilon());
    ASSERT_DOUBLE_EQ(-1., pts1.back());

    auto pts2 = detail::chebyshev_points<double>(0);
    ASSERT_EQ(pts2.size(), 0);
}

TEST(ChebfunPoly, FirstFewPolys) {
    size_t N = 7;
    auto pts = detail::chebyshev_points<double>(N);

    std::vector<double> values(N);

    // test x^2
    for (size_t i = 0; i < N; ++i)
        values[i] = pts[i]*pts[i];

    auto res = chebpoly(values);

    // Chebyshev coefficients should be [0.5, 0, 0.5, ....]
    for (size_t i = 0; i < N; ++i) {
        if (i == 0 || i == 2)
            ASSERT_LT(std::abs(res[i] - 0.5), std::numeric_limits<double>::epsilon());
        else
            ASSERT_LT(std::abs(res[i]), std::numeric_limits<double>::epsilon());
    }

    // test x^3
    for (size_t i = 0; i < N; ++i)
        values[i] = pts[i]*pts[i]*pts[i];

    res = chebpoly(values);

    // Chebyshev coefficients should be [0, 0.75, 0, 0.25, ....]
    for (size_t i = 0; i < N; ++i) {
        if (i == 1)
            ASSERT_LT(std::abs(res[i] - 0.75), std::numeric_limits<double>::epsilon());
        else if (i == 3)
            ASSERT_LT(std::abs(res[i] - 0.25), std::numeric_limits<double>::epsilon());
        else
            ASSERT_LT(std::abs(res[i]), std::numeric_limits<double>::epsilon());
    }
}

/// Test a sample function call.
///
#if 0
TEST(add, test)
{
 //ASSERT_EQ(3, add(1, 2));
}
#endif
