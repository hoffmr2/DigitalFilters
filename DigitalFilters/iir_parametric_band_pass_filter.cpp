#include "iir_parametric_band_pass_filter.h"



IIRParametricBandPassFilter::IIRParametricBandPassFilter(double cutoff_frequency, double sample_rate, double gain, double q_factor, bypassState bypass_state) : DigitalFilter(cutoff_frequency,sample_rate,bypass_state),
																																								q_factor_(q_factor),
																																								gain_flag_(true)
{
	SetAngularCutoffFrequency();
	CalculateBandwidth();
	SetGainDb(gain);
	InitMemory();
	InitAcoefficients();
	InitBcoefficients();
}

IIRParametricBandPassFilter::~IIRParametricBandPassFilter()
{
}

void IIRParametricBandPassFilter::CalculateBandwidth()
{
	assert(q_factor_ >= 0);
	bandwidth_ = q_factor_*angular_cutoff_frequency_;
}

void IIRParametricBandPassFilter::ChangeCutoffFrequency(double cutoff_frequency)
{
	assert(cutoff_frequency > 0);
	cutoff_frequency_ = cutoff_frequency;

	SetAngularCutoffFrequency();
	CalculateBandwidth();
	InitAcoefficients();
	InitBcoefficients();
}

void IIRParametricBandPassFilter::ChangeQFactor(double q_factor)
{
	assert(q_factor >= 0);
	q_factor_ = q_factor;

	SetAngularCutoffFrequency();
	CalculateBandwidth();
	InitAcoefficients();
	InitBcoefficients();
}

float IIRParametricBandPassFilter::FilterOutputLeft(float sample)
{
	if (bypass_state_ == off)
		return FilterOutput(sample, memory_left_);
	else
		return sample;

}

void IIRParametricBandPassFilter::CalculateCorrectedBCoefficients(double& b0, double& b1, double& b2) const
{
	b0 = b_coefficients_[0] +gain_*b_coefficients_[1];
	b1 = b_coefficients_[2];
	b2 = b_coefficients_[0] - gain_*b_coefficients_[1];
}

float IIRParametricBandPassFilter::FilterOutputRight(float sample)
{
	if (bypass_state_ == off)
		return FilterOutput(sample, memory_right_);
	else
		return sample;
}

float IIRParametricBandPassFilter::FilterOutput(float sample, float* memory) const
{
	double ans;
	double b0,b1,b2;	
	CalculateCorrectedBCoefficients(b0, b1, b2);


	if (gain_flag_ == true)
	{
		auto tmp = (sample - a_coefficients_[1] * memory[0] - a_coefficients_[2] * memory[1]);
		ans = tmp*b0 + b1*memory[0] + b2*memory[1];
		memory[1] = memory[0];
		memory[0] = float (tmp);
	}
	else
	{
		auto tmp = (sample - b1 * memory[0] / b0 - b2 * memory[1] / b0);
		ans = tmp / b0 + a_coefficients_[1] * memory[0] / b0 + a_coefficients_[2] * memory[1] / b0;
		memory[1] = memory[0];
		memory[0] = float (tmp);
	}

	return float(ans);
}

double IIRParametricBandPassFilter::Spectrum(double frequency)
{
	if (bypass_state_ == off)
	{
		double b0, b1, b2;
		CalculateCorrectedBCoefficients(b0, b1, b2);

		if (gain_flag_ == true)
			return sqrt((pow((b0 + b1*cos(2 * PI*frequency) + b2*cos(4 * PI*frequency)), 2) + pow(b1*sin(2 * PI*frequency) + b2*sin(4 * PI*frequency), 2)) /
			(pow((1 + a_coefficients_[1] * cos(2 * PI*frequency) + a_coefficients_[2] * cos(4 * PI*frequency)), 2) + pow(a_coefficients_[1] * sin(2 * PI*frequency) + a_coefficients_[2] * sin(4 * PI*frequency), 2)));
		else
			return sqrt((pow((1 + a_coefficients_[1] * cos(2 * PI*frequency) + a_coefficients_[2] * cos(4 * PI*frequency)), 2) + pow(a_coefficients_[1] * sin(2 * PI*frequency) + a_coefficients_[2] * sin(4 * PI*frequency), 2)) /
			(pow((b0 + b1*cos(2 * PI*frequency) + b2*cos(4 * PI*frequency)), 2) + pow(b1*sin(2 * PI*frequency) + b2*sin(4 * PI*frequency), 2)));
	}
	else
		return 1;
}

void IIRParametricBandPassFilter::InitAcoefficients()
{
	
	if (a_coefficients_ == nullptr)
		a_coefficients_ = new double[COEFFICIENTS_NUMBER];

	auto denominator = 4 * sample_rate_*sample_rate_ + angular_cutoff_frequency_*angular_cutoff_frequency_ + sample_rate_*(2 * bandwidth_);
	a_coefficients_[0] = 1;
	a_coefficients_[1] = (angular_cutoff_frequency_*(2 * angular_cutoff_frequency_) / denominator) - sample_rate_*((8 * sample_rate_) / denominator);
	a_coefficients_[2] = ((4 * sample_rate_) / denominator)*sample_rate_ + angular_cutoff_frequency_*(angular_cutoff_frequency_ / denominator) - bandwidth_*((2 * sample_rate_) / denominator);
}

void IIRParametricBandPassFilter::InitBcoefficients()
{
	if (b_coefficients_ == nullptr)
		b_coefficients_ = new double[COEFFICIENTS_NUMBER];

	auto denominator = 4 * sample_rate_*sample_rate_ + angular_cutoff_frequency_*angular_cutoff_frequency_ + 2 * bandwidth_*sample_rate_;
	b_coefficients_[0] = ((4 * sample_rate_*sample_rate_) / denominator) + ((angular_cutoff_frequency_*angular_cutoff_frequency_) / denominator);
	b_coefficients_[1] = 2 * bandwidth_*(sample_rate_ / denominator);
	b_coefficients_[2] = ((2 * angular_cutoff_frequency_*angular_cutoff_frequency_) / denominator) - ((8 * sample_rate_*sample_rate_) / denominator);
}

void IIRParametricBandPassFilter::InitMemory()
{
	if (memory_left_ == nullptr)
		memory_left_ = new float[MEMORY_SIZE];
	if (memory_right_ == nullptr)
		memory_right_ = new float[MEMORY_SIZE];
	for(auto i=0;i<MEMORY_SIZE;++i)
	{
		memory_left_[i] = 0;
		memory_right_[i] = 0;
	}
}

void IIRParametricBandPassFilter::SetGainDb(double gain_db)
{
	gain_flag_ = (gain_db <= 0) ? true : false;
	auto exp = (gain_db) / 20.0;
	gain_ = std::pow(10, (-std::abs(exp)));
}
