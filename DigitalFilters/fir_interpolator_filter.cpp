#include "fir_interpolator_filter.h"






FIRInterpolatorFilter::FIRInterpolatorFilter(int filter_size, int interpolation_factor) : DigitalFilter(1,1), filter_size_(filter_size), interpolation_factor_(interpolation_factor)
{
	InitMemory();
	InitBcoefficients();
}

FIRInterpolatorFilter::~FIRInterpolatorFilter()
{
	if (b_coefficients_ != nullptr)
		delete[] b_coefficients_;
	if (memory_left_ != nullptr)
		delete[] memory_left_;
}

float FIRInterpolatorFilter::FilterOutput(float* samples, int start, int samples_vector_size) const
{
	float ans =0;
	for (int i = 0; i < filter_size_; ++i)
			ans += float(memory_left_[i] * b_coefficients_[i]);
	ans += samples[start];
	ShiftMemory(samples[start]);
	for (int i = start +1; i <start + filter_size_ && i < samples_vector_size; ++i)
		ans += float(samples[i] * b_coefficients_[i-start]);

	return ans;
}

void FIRInterpolatorFilter::InitBcoefficients()
{
	b_coefficients_ = new double[2 * filter_size_ + 1];

	for (auto i=0;i<2*filter_size_+1;++i)
	{
		if (i - filter_size_ != 0)
			b_coefficients_[i] = sin(PI*(i - filter_size_) / interpolation_factor_) / (PI*(i - filter_size_) / interpolation_factor_);
		else
			b_coefficients_[i] = 1;
		auto hanning_window = (0.54 - 0.46* cos(2 * PI*(i-filter_size_) / (2 * filter_size_)));
		b_coefficients_[i] *= hanning_window;
	}
}

void FIRInterpolatorFilter::InitMemory()
{
	if (memory_left_ == nullptr)
		memory_left_ = new float[filter_size_];
	for (auto i = 0; i < filter_size_; ++i)
		memory_left_[i] = 0;
}

void FIRInterpolatorFilter::ShiftMemory(float sample) const
{
	float tmp;
	for (auto i = filter_size_-2; i >=0 ; --i)
	{
		tmp = memory_left_[i];
		memory_left_[i] = memory_left_[i + 1];
	}
	memory_left_[filter_size_ - 1] = sample;
}
