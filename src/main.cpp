#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include "AudioFile.h"
#include <vector>
#include <numeric>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
using namespace matplot;

struct Wave {
    std::vector<double> x;
    std::vector<double> y;
};

const double PI = 3.14;

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

void generate_sine_wave(double frequency) {
    Wave sin;
    for (int i = 0; i < 2*PI*100; i++){
        sin.x.push_back(static_cast<double>(i) / (100*PI));
        sin.y.push_back(std::sin(PI * frequency * sin.x[i]));
    }
    title("Sine wave");
    plot(sin.x, sin.y);
    show();
}

void generate_cosine_wave(double frequency) {
    Wave cos;
    for (int i = 0; i < 628; i++){
        cos.x.push_back(static_cast<double>(i) / (100*PI));
        cos.y.push_back(std::cos(3.14 * frequency * cos.x[i]));
    }
    title("Cosine wave");
    plot(cos.x, cos.y);
    show();
}

void generate_sawtooth_wave(double frequency) {

}

void generate_square_wave(double frequency) {


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

    m.def("generate_sine_wave", &generate_sine_wave, R"pbdoc(
        
    generate and display sine wave with given frequency    

    )pbdoc");

    m.def("generate_cosine_wave", &generate_cosine_wave, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");
    

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
