#include "fir_polyphase_decimator_filter.h"





FirPolyphaseDecimatorFilter::FirPolyphaseDecimatorFilter(double sample_rate, int decimation_factor)
	:FirLowPassFilter(sample_rate,sample_rate/(2*decimation_factor),sample_rate/decimation_factor,60), decimation_factor_(decimation_factor),polyphase_filter_size_(decimation_factor),
	sample_counter_left_(-1),sample_counter_right_(-1)
{
	InitDecimatorFilter();
}

FirPolyphaseDecimatorFilter::~FirPolyphaseDecimatorFilter()
{
}

void FirPolyphaseDecimatorFilter::CreatePolyphaseFilter(double* tmp_b_coefficients, int i)
{
	for(int j=0;j<polyphase_filter_size_;++j)
		if (i + j*decimation_factor_ < filter_size_)
			tmp_b_coefficients[j] = b_coefficients_[i + j*decimation_factor_];
		else
			tmp_b_coefficients[j] = 0;

	polyphase_filters_.push_back(FirFilter(sample_rate_ / decimation_factor_, 1, tmp_b_coefficients, polyphase_filter_size_));
}

void FirPolyphaseDecimatorFilter::InitPolyphaseFilters()
{
	double* tmp_b_coefficients = new double[polyphase_filter_size_];
	for (int i = 0; i < number_of_polyphase_filters; ++i)
	{
		
		CreatePolyphaseFilter(tmp_b_coefficients, i);
	}
	delete [] tmp_b_coefficients;
}

void FirPolyphaseDecimatorFilter::InitDecimatorFilter()
{

	assert(filter_size_ != 0);
	assert(polyphase_filter_size_ != 0);

	number_of_polyphase_filters = int(ceil(filter_size_ / decimation_factor_));
	InitPolyphaseFilters();
}

float FirPolyphaseDecimatorFilter::FilterOutputLeft(float sample)
{
	++sample_counter_left_;
	sample_counter_left_ %= decimation_factor_;
	return FilterOutput(sample, sample_counter_left_);
	
}

float FirPolyphaseDecimatorFilter::FilterOutputRight(float sample)
{
	++sample_counter_right_;
	sample_counter_right_ %= decimation_factor_;
	return FilterOutput(sample, sample_counter_right_);
}

float FirPolyphaseDecimatorFilter::FilterOutput(float sample, int& sample_counter)
{
	if (sample_counter == 0)
		return polyphase_filters_[sample_counter].FilterOutputLeft(sample);
	else
		return polyphase_filters_[polyphase_filter_size_ - sample_counter - 1].FilterOutputLeft(sample);
}
