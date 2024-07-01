
#include "idaklu.hpp" 

#include "idaklu_casadi.hpp" 
#include "idaklu_python.hpp" 

#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <vector>

Function generate_function(const std::string& data) {
  return Function::deserialize(data);
}

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<np_array>);

PYBIND11_MODULE(idaklu, m)
{
  m.doc() = "sundials solvers"; // optional module docstring

  py::bind_vector<std::vector<np_array>>(m, "VectorNdArray");

  m.def("solve_python", &solve_python, "The solve function for python evaluators", 
        py::arg("t"), py::arg("y0"),
        py::arg("yp0"), py::arg("res"), py::arg("jac"), py::arg("sens"), 
        py::arg("get_jac_data"),
        py::arg("get_jac_row_vals"), py::arg("get_jac_col_ptr"), py::arg("nnz"),
        py::arg("events"), py::arg("number_of_events"), py::arg("use_jacobian"),
        py::arg("rhs_alg_id"), py::arg("atol"), py::arg("rtol"), py::arg("inputs"),
        py::arg("number_of_sensitivity_parameters"),
        py::return_value_policy::take_ownership);

  m.def("solve_casadi", &solve_casadi, "The solve function for casadi evaluators", 
        py::arg("t"), py::arg("y0"), py::arg("yp0"), 
        py::arg("rhs_alg"), 
        py::arg("jac_times_cjmass"), 
        py::arg("jac_times_cjmass_colptrs"), 
        py::arg("jac_times_cjmass_rowvals"), 
        py::arg("jac_times_cjmass_nnz"), 
        py::arg("jac_action"), 
        py::arg("mass_action"), 
        py::arg("sens"), 
        py::arg("events"), py::arg("number_of_events"), 
        py::arg("use_jacobian"),
        py::arg("rhs_alg_id"),
        py::arg("atol"), py::arg("rtol"), py::arg("inputs"),
        py::arg("number_of_sensitivity_parameters"),
        py::return_value_policy::take_ownership);

  m.def("generate_function", &generate_function, "Generate a casadi function", 
        py::arg("string"),
        py::return_value_policy::take_ownership);

  py::class_<Function>(m, "Function");

  py::class_<Solution>(m, "solution")
      .def_readwrite("t", &Solution::t)
      .def_readwrite("y", &Solution::y)
      .def_readwrite("yS", &Solution::yS)
      .def_readwrite("flag", &Solution::flag);
}


