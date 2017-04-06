//
// Created by Anton Loukianov on 3/30/17.
//

#include "chebfun_module.h"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

PYBIND11_PLUGIN(chebfunpy)
{
    py::module m("chebfunpy", "a toy chebfun implementation in C++");

    py::class_<chebfun::basic_chebfun<double>> basic_chebfun(m, "BasicChebfun");

    m.def("chebpoly", &chebfun::chebpoly<double>);
    m.def("chebpts", &chebfun::detail::chebyshev_points<double>);
    return m.ptr();
}


