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
    std::vector<double> y_imag;
    std::vector<std::complex<double>> x_complex;
};

const double PI = 3.14;
const int BUFF = 32767;
typedef std::complex<double> Complex;

void generate_sine_wave(double frequency) {
    Wave sin;
    for (int i = 0; i < 2*PI*100; i++){
        sin.x.push_back(static_cast<double>(i) / (100*PI));
        sin.y.push_back(std::sin(PI * frequency * sin.x[i]));
    }
    plot(sin.x, sin.y);
    show();
}

void generate_cosine_wave(double frequency) {
    Wave cos;
    for (int i = 0; i < 2*PI*100; i++){
        cos.x.push_back(static_cast<double>(i) / (100*PI));
        cos.y.push_back(std::cos(PI * frequency * cos.x[i]));
    }
    plot(cos.x, cos.y);
    show();
}

void visualize_audio(std::string audio_path) {
    AudioFile<double> audio_file;
    Wave audio_wave;
    bool loaded = audio_file.load(audio_path);
    if (!loaded) {
        std::cerr << "Audio is not loaded!" << std::endl;
        return;
    }
    else {
        int num_of_samples = 500;
        for (int i = 0; i < num_of_samples; i++) {
            audio_wave.x.push_back(static_cast<double>(i));
            audio_wave.y.push_back(static_cast<double>(audio_file.samples[0][i]));
        }
    }
    plot(audio_wave.x, audio_wave.y);
    show();
}

void generate_sawtooth_wave(double frequency, int length) {
    Wave saw; 
    double amplitude = 1.0;

    for (int i = 0; i < length; ++i) {
        double t = static_cast<double>(i) / length;
        saw.y.push_back(amplitude * (2.0 * (t * frequency - std::floor(0.5 + t * frequency))));
    }

    std::iota(saw.x.begin(), saw.x.end(), 0);

    plot(saw.x, saw.y);
    title("Sawtooth Wave");
    xlabel("Sample Index");
    ylabel("Amplitude");

    show();   
}

void generate_square_wave(double frequency, int length, int prog) {
    Wave square;

    int period = length / (2 * frequency);

    for (size_t i = 0; i < length; ++i) {

        if ((i / period) % 2 == 0) {
            square.y.push_back(BUFF);
        }
        else {
            square.y.push_back(-BUFF);
        }
    }
    std::iota(square.x.begin(), square.x.end(), 0); 

    plot(square.x, square.y);
    title("Square Wave");
    xlabel("Sample Index");
    ylabel("Amplitude");
    show();
}

void threshold_signal(std::string audio_path,double threshhold) {
	AudioFile<double> audio_file;
	Wave audio_wave;
	bool loaded = audio_file.load(audio_path);

	int num_of_samples = 500;
		for (int i = 0; i < num_of_samples; i++) {
			audio_wave.x.push_back(static_cast<double>(i));
			audio_wave.y.push_back(static_cast<double>(audio_file.samples[0][i]));
		}
	
	for (int i = 0; i < 500; i++) {
		if (audio_wave.y.at(i) > threshhold) {
			audio_wave.y.at(i) = 1;
		}
		else {
			audio_wave.y.at(i) = 0;
		}
	}
	plot(audio_wave.x, audio_wave.y);
	title("threshold signal");
	show();
}


void compute_and_plot_dft(double frequency, double amplitude, double sampleRate, int numSamples) {
    Wave begin_wave;
    Wave dft_wave;
    Wave idft_wave;

    for (int i = 0; i < numSamples; ++i) {
        double t = i / sampleRate;
        double sawtooth_value = amplitude * (2 * (t * frequency - floor(t * frequency + 0.5)));
        begin_wave.x.push_back(t);
        begin_wave.y.push_back(sawtooth_value);
    }
    plot(begin_wave.x, begin_wave.y);
    show();

    int N = begin_wave.y.size();

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = -2.0 * PI * k * n / N;
            std::complex<double> w(std::cos(angle), std::sin(angle));
            sum += begin_wave.y[n] * w;
        }
        double t = k / sampleRate;
        dft_wave.x.push_back(t);
        dft_wave.x_complex.push_back(sum);
        dft_wave.y.push_back(std::abs(sum));
    }
    plot(dft_wave.x, dft_wave.y);
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

    plot(dft_wave.x, idft_wave.y);
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

    m.def("generate_sine_wave", &generate_sine_wave, R"pbdoc(
        
    generate and display sine wave with given frequency    

    )pbdoc");

    m.def("generate_cosine_wave", &generate_cosine_wave, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");


    m.def("threshold_signal", &threshold_signal, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");

    m.def("visualize_audio", &visualize_audio, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");

    m.def("generate_sawtooth_wave", &generate_sawtooth_wave, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");

    m.def("generate_square_wave", &generate_square_wave, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");

    m.def("compute_and_plot_dft", &compute_and_plot_dft, R"pbdoc(
        
    generate and display cosine wave with given frequency    

    )pbdoc");
    

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
