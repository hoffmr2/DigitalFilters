#include "fir_polyphase_decimator_filter.h"





FirPolyphaseDecimatorFilter::FirPolyphaseDecimatorFilter(double sample_rate, int decimation_factor)
	:FirLowPassFilter(sample_rate, sample_rate / (3 * decimation_factor), sample_rate / decimation_factor, 80), decimation_factor_(decimation_factor), number_of_polyphase_filters(decimation_factor),
	sample_counter_left_(-1), sample_counter_right_(-1)
{
	InitDecimatorFilter();
}

FirPolyphaseDecimatorFilter::~FirPolyphaseDecimatorFilter()
{
}

void FirPolyphaseDecimatorFilter::CreatePolyphaseFilter(double* tmp_b_coefficients, int i)
{
	for (int j = 0; j < polyphase_filter_size_; ++j)
	{
		if (i + j*decimation_factor_ < filter_size_)
			tmp_b_coefficients[j] = b_coefficients_[i + j*decimation_factor_];
		else
			tmp_b_coefficients[j] = 0;
	}
	polyphase_filters_.push_back(new FirFilter(sample_rate_ / decimation_factor_, 1, tmp_b_coefficients, polyphase_filter_size_));
}

void FirPolyphaseDecimatorFilter::InitPolyphaseFilters()
{
	double* tmp_b_coefficients = new double[polyphase_filter_size_];
	for (int i = 0; i < number_of_polyphase_filters; ++i)
	{

		CreatePolyphaseFilter(tmp_b_coefficients, i);
	}
	delete[] tmp_b_coefficients;
}

void FirPolyphaseDecimatorFilter::InitFilterOutpuTables()
{
	output_left_ = new float[number_of_polyphase_filters];
	output_right_ = new float[number_of_polyphase_filters];

	for (int i = 0; i<number_of_polyphase_filters; ++i)
	{
		output_left_[i] = 0;
		output_right_[i] = 0;
	}
}

void FirPolyphaseDecimatorFilter::InitDecimatorFilter()
{

	assert(filter_size_ != 0);


	polyphase_filter_size_ = int(ceil(double(filter_size_) / decimation_factor_));

	InitFilterOutpuTables();


	InitPolyphaseFilters();
}

float FirPolyphaseDecimatorFilter::FilterOutputLeft(float sample)
{
	++sample_counter_left_;
	sample_counter_left_ %= decimation_factor_;
	if (sample_counter_left_ == 0)
		output_right_[sample_counter_left_] = polyphase_filters_[sample_counter_left_]->FilterOutputLeft(sample);
	else
	{
		int index = number_of_polyphase_filters - sample_counter_left_;
		output_left_[index] = polyphase_filters_[index]->FilterOutputLeft(sample);
	}

	float ans = 0;
	for (int i = 0; i < number_of_polyphase_filters; ++i)
		ans += output_left_[i];

	return ans;

}

float FirPolyphaseDecimatorFilter::FilterOutputRight(float sample)
{
	++sample_counter_right_;
	sample_counter_right_ %= decimation_factor_;
	if (sample_counter_right_ == 0)
		output_right_[sample_counter_right_] = polyphase_filters_[sample_counter_right_]->FilterOutputRight(sample);
	else
	{
		int index = number_of_polyphase_filters - sample_counter_right_ - 1;
		output_right_[index] = polyphase_filters_[index]->FilterOutputRight(sample);
	}

	float ans = 0;
	for (int i = 0; i < number_of_polyphase_filters; ++i)
		ans += output_right_[i];

	return ans;
}

float FirPolyphaseDecimatorFilter::FilterOutput(float sample, int& sample_counter)
{
	if (sample_counter == 0)
		return polyphase_filters_[sample_counter]->FilterOutputLeft(sample);
	else
		return polyphase_filters_[polyphase_filter_size_ - sample_counter - 1]->FilterOutputLeft(sample);
}


