#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int a, int b) { return a + b; }

#define PYBIND11_MODULE_(name, m) PYBIND11_MODULE(name, m)

PYBIND11_MODULE_(PYMPSYM, m) {
  m.doc() = DESCRIPTION;

  m.def("add", &add, "Add two numbers"); // TODO
}
