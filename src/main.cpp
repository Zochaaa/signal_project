#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <AudioFile.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <string>
#include <pybind11/complex.h>
#include <complex>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
using namespace matplot;

struct Wave {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<std::complex<double>> x_complex;
    double frequency;
    std::string audio_path;
    int length;
};

const double PI = 3.14;
const int BUFF = 32767;
typedef std::complex<double> Complex;

Wave generate_sine_wave(double frequency) {
    Wave sin;
    sin.length = 2 * PI * 100;
    sin.frequency = frequency;
    for (int i = 0; i < 2*PI*100; i++){
        sin.x.push_back(static_cast<double>(i) / (100*PI));
        sin.y.push_back(std::sin(PI * frequency * sin.x[i]));
    }
    plot(sin.x, sin.y);
    show();
    return sin;
}

Wave generate_cosine_wave(double frequency) {
    Wave cos;
    cos.length = 2 * PI * 100;
    cos.frequency = frequency;
    for (int i = 0; i < 2*PI*100; i++){
        cos.x.push_back(static_cast<double>(i) / (100*PI));
        cos.y.push_back(std::cos(PI * frequency * cos.x[i]));
    }
    plot(cos.x, cos.y);
    show();
    return cos;
}

Wave visualize_audio(std::string audio_path) {
    AudioFile<double> audio_file;
    Wave audio_wave;
    
    bool loaded = audio_file.load(audio_path);
    if (!loaded) {
        std::cerr << "Audio is not loaded!" << std::endl;
        visualize_audio(audio_path);
    }
    else {
        audio_wave.audio_path = audio_path;
        double time = audio_file.getLengthInSeconds();
        int sample_rate = audio_file.getSampleRate();
        audio_wave.length = time * sample_rate;
        for (int i = 0; i < time*sample_rate; i++) {
            audio_wave.y.push_back(static_cast<double>(audio_file.samples[0][i]));
        }
        for (int i = 0; i <= time; i++) {
            audio_wave.x.push_back(i);
        }
    }
    plot(audio_wave.x, audio_wave.y);
    show();
    return audio_wave;
}

Wave generate_sawtooth_wave(double frequency, int length) {
    Wave saw; 
    saw.length = length;
    saw.frequency = frequency;
    double amplitude = 1.0;
    for (int i = 0; i < length; ++i) {
        double t = static_cast<double>(i) / length;
        saw.y.push_back(amplitude * (2.0 * (t * frequency - std::floor(0.5 + t * frequency))));
        saw.x.push_back(t);
    }
    //std::iota(saw.x.begin(), saw.x.end(a), 0);

    plot(saw.x, saw.y);
    show();   
    return saw;
}

Wave generate_square_wave(double frequency, int length) {
    Wave square;
    square.length = length;
    square.frequency = frequency;
    int period = length / (2 * frequency);

    for (size_t i = 0; i < length; ++i) {

        if ((i / period) % 2 == 0) {
            square.y.push_back(BUFF);
           
        }
        else {
            square.y.push_back(-BUFF);
        }
        square.x.push_back(i);
    }
    //std::iota(square.x.begin(), square.x.end(), 0); 

    plot(square.x, square.y);
    show();
    return square;
}

void threshold_signal(Wave &begin_wave,double threshhold) {
    Wave end_wave;
	
	for (int i = 0; i < begin_wave.length; i++) {
		if (end_wave.y.at(i) > threshhold) {
			end_wave.y.at(i) = 1.0;
		}
		else {
			end_wave.y.at(i) = 0.0;
		}
	}
	plot(end_wave.x, end_wave.y);
	show();
}


void dft_idft(double frequency, double amplitude, double sample_rate, int num_samples) {
    Wave begin_wave;
    Wave dft_wave;
    Wave idft_wave;

    for (int i = 0; i < num_samples; ++i) {
        double t = i / sample_rate;
        double sawtooth_value = amplitude * (2 * (t * frequency - floor(t * frequency + 0.5)));
        begin_wave.x.push_back(t);
        begin_wave.y.push_back(sawtooth_value);
    }
    auto ax1 = subplot(2, 2, 0);
    plot(ax1, begin_wave.x, begin_wave.y);\
    title(ax1, "Sawtooth");
    show();

    int N = begin_wave.y.size();

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * PI * k * n / N;
            std::complex<double> w(std::cos(angle), std::sin(angle));
            sum += begin_wave.y[n] * w;
        }
        double t = k / sample_rate;
        dft_wave.x.push_back(t);
        dft_wave.x_complex.push_back(sum);
        dft_wave.y.push_back(std::abs(sum));
    }
    auto ax2 = subplot(2, 2, 1);
    plot(ax2, dft_wave.x, dft_wave.y);
    title(ax2, "DFT");
    show();

    N = dft_wave.x_complex.size();
    idft_wave.y.resize(N);

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = 2.0 * PI * k * n / N;
            std::complex<double> w(std::cos(angle), std::sin(angle));
            sum += dft_wave.x_complex[n] * w;
        }
		idft_wave.y[k] = std::real(sum) / N;
    }
    auto ax3 = subplot(2, 2, 2);
    plot(ax3, dft_wave.x, idft_wave.y);
    title(ax3, "IDFT");
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

         
    )pbdoc";

    py::class_<Wave>(m, "Wave")
        .def(py::init<>())
        .def_readwrite("x", &Wave::x)
        .def_readwrite("y", &Wave::y)
        .def_readwrite("y_imag", &Wave::y_imag)
        .def_readwrite("x_complex", &Wave::x_complex)
        .def_readwrite("frequency", &Wave::frequency)
        .def_readwrite("audio_path", &Wave::audio_path)
        .def_readwrite("length", &Wave::length);

    m.def("generate_sine_wave", &generate_sine_wave, R"pbdoc(
        
    generate and display sine wave with given frequency    

    )pbdoc");

    m.def("generate_cosine_wave", &generate_cosine_wave, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");


    m.def("threshold_signal", &threshold_signal, R"pbdoc(
        
      

    )pbdoc");

    m.def("visualize_audio", &visualize_audio, R"pbdoc(
        
    generate  wave with audio  

    )pbdoc");

    m.def("generate_sawtooth_wave", &generate_sawtooth_wave, R"pbdoc(
        
    generate and display sawtooth wave with given frequency    

    )pbdoc");

    m.def("generate_square_wave", &generate_square_wave, R"pbdoc(
        
    generate and display square wave with given frequency    

    )pbdoc");

    m.def("dft_idft", &dft_idft, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");
    

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
