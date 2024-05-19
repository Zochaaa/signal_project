#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <vector>
#include <numeric>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
using namespace matplot;

int add(int i, int j) {
    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {2, 3, 5, 7, 11};

    // Plot data
    plot(x, y);
    title("Simple Plot");
    xlabel("X-axis");
    ylabel("Y-axis");

    // Show the plot
    show();
    return i + j;
}

int minus(int i, int j) {
    return i - j;
}
void plot(){
    int length = 1000;
    double frequency = 5.0;
    double amplitude = 1.0;
     // Generate the sawtooth signal
    std::vector<double> signal(length);
    for (int i = 0; i < length; ++i) {
        double t = static_cast<double>(i) / length;
        signal[i] = amplitude * (2.0 * (t * frequency - std::floor(0.5 + t * frequency)));
    }

    // Create x values
    std::vector<double> x(signal.size());
    std::iota(x.begin(), x.end(), 0);

    // Plot the signal
    plot(x, signal);
    title("Sawtooth Signal");
    xlabel("Sample Index");
    ylabel("Amplitude");

    // Show the plot
    show();
} 

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           minus
           subtract
           plot
    )pbdoc";


    m.def("plot", &plot, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    
    m.def("minus", &minus, R"pbdoc(
        Subtruct two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
