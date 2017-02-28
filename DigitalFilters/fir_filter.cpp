#include "fir_filter.h"



FirFilter::FirFilter(double sample_rate, double pass_band_frequency)
	: DigitalFilter(pass_band_frequency, sample_rate), filter_size_(0)

{


}

FirFilter::FirFilter(double sample_rate, double pass_band_frequency, double* b_coefficients, int filter_size)
	: DigitalFilter(pass_band_frequency, sample_rate), filter_size_(filter_size), tmp_b_coeffcients_(b_coefficients)
{
	InitMemory();
	InitBcoefficients();
}

FirFilter::~FirFilter()
{
}

float FirFilter::FilterOutputLeft(float sample)
{
	return FilterOutput(sample, memory_left_);
}

float FirFilter::FilterOutputRight(float sample)
{
	return FilterOutput(sample, memory_right_);
}

double FirFilter::Spectrum(double frequency)
{
	assert(b_coefficients_ != nullptr);
	auto real = 0.0;
	auto imag = 0.0;

	for (auto i = 0; i<filter_size_; ++i)
	{
		real += cos(2 * PI*frequency / sample_rate_)*b_coefficients_[i];
		imag += sin(2 * PI*frequency / sample_rate_)*b_coefficients_[i];
	}

	return sqrt(real*real + imag*imag);
}

void FirFilter::ChangeCutoffFrequency(double newFpass)
{
}

void FirFilter::InitBcoefficients()
{
	assert(filter_size_ != 0);
	if (b_coefficients_ != nullptr)
		delete b_coefficients_;
	b_coefficients_ = new double[filter_size_];

	for (int i = 0; i < filter_size_; ++i)
		b_coefficients_[i] = tmp_b_coeffcients_[i];

}

void FirFilter::InitAcoefficients()
{
	a_coefficients_ = nullptr;
}




void FirFilter::ShiftMemory(float* memory) const
{
	auto tmp = memory[0];
	for (auto i = 1; i<filter_size_; ++i)
	{
		memory[i] = memory[i - 1];
	}
}

void FirFilter::ClearMemory() const
{
	for (auto i = 0; i<filter_size_; ++i)
	{
		memory_left_[i] = 0;
		memory_right_[i] = 0;
	}
}

void FirFilter::InitMemory()
{
	assert(filter_size_ != 0);
	if (memory_left_ == nullptr &&
		memory_right_ == nullptr)
	{
		memory_left_ = new float[filter_size_];
		memory_right_ = new float[filter_size_];
	}
	ClearMemory();
}

float FirFilter::FilterOutput(float sample, float* memory) const
{
	assert(memory != nullptr);
	assert(b_coefficients_ != nullptr);

	memory[0] = sample;
	double output = 0;

	for (auto i = 0; i < filter_size_; ++i)
		output += memory[i] * b_coefficients_[i];

	ShiftMemory(memory);

	return output;
}
