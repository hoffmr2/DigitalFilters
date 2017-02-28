#pragma once
#include "fir_low_pass_filter.h"
#include <vector>


class FirPolyphaseDecimatorFilter : FirLowPassFilter
{
public:
	FirPolyphaseDecimatorFilter(double sample_rate, int decimation_factor);
	~FirPolyphaseDecimatorFilter();
	void CreatePolyphaseFilter(double* tmp_b_coefficients, int i);
	void InitPolyphaseFilters();
	void InitFilterOutpuTables();
	void InitDecimatorFilter();

	virtual float FilterOutputLeft(float sample) override;
	virtual float FilterOutputRight(float sample) override;

private:
	float FilterOutput(float sample, int& sample_counter);
	float* output_left_;
	float* output_right_;
	const int decimation_factor_;
	std::vector<FirFilter*> polyphase_filters_;
	int number_of_polyphase_filters;
	int polyphase_filter_size_;
	int sample_counter_left_;
	int sample_counter_right_;

};

