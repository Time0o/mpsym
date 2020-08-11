#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int a, int b) { return a + b; }

PYBIND11_MODULE(mpsym, m) {
  m.doc() = "MPSoC symmetry reduction package";

  m.def("add", &add, "Add two numbers"); // TODO
}
