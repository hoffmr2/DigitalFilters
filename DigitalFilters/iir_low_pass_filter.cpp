#include "iir_low_pass_filter.h"


IIRLowPassFilter::IIRLowPassFilter(double cutoff_frequency, double sample_rate, double absorbtion_db, bypassState bypass_state)
	:DigitalFilter(cutoff_frequency,sample_rate,bypass_state)
{
	assert(absorbtion_db > 0);
	SetAbsorbtionFactor(absorbtion_db);
	InitMemory();
	InitAcoefficients();
	InitBcoefficients();
}

IIRLowPassFilter::~IIRLowPassFilter()
{
}

void IIRLowPassFilter::ChangeCutoffFrequency(double newFpass)
{
	assert(newFpass > 0);
	cutoff_frequency_ = newFpass;
	SetAngularCutoffFrequency();
	InitMemory();
	InitAcoefficients();
	InitBcoefficients();

}

double IIRLowPassFilter::GetCoefficientDenominator() const
{
	return 4 * absorbtion_factor_ * absorbtion_factor_ * sample_rate_ * sample_rate_ +
		2 * absorbtion_factor_ * angular_cutoff_frequency_ * sample_rate_ * sqrt(2) +
		angular_cutoff_frequency_ * angular_cutoff_frequency_;
}

double IIRLowPassFilter::GetA1CoefficientValue() const
{
	auto denominator = GetCoefficientDenominator();
	return (2 * angular_cutoff_frequency_ * angular_cutoff_frequency_ - 
		8 * absorbtion_factor_ * absorbtion_factor_ * sample_rate_ * sample_rate_) / denominator;
}

int IIRLowPassFilter::GetA2Coefficient() const
{
	auto denominator = GetCoefficientDenominator();
	return (angular_cutoff_frequency_ * angular_cutoff_frequency_ + 
		4 * absorbtion_factor_ * absorbtion_factor_ * sample_rate_ * sample_rate_ - 
			2 * absorbtion_factor_ * angular_cutoff_frequency_ * sample_rate_ * sqrt(2)) / denominator;
}

float IIRLowPassFilter::FilterOutputLeft(float sample)
{
	if (bypass_state_ == off)
		return FilterOutput(memory_right_, sample);
	else
		return sample;
}

float IIRLowPassFilter::FilterOutputRight(float sample)
{
	if (bypass_state_ == off)
		return FilterOutput(memory_left_, sample);
	else
		return sample;
}

double IIRLowPassFilter::Spectrum(double frequency)
{
	return sqrt((pow((b_coefficients_[0] + b_coefficients_[1]*cos(2 * PI*frequency) + b_coefficients_[2]*cos(4 * PI*frequency)), 2)
		+ pow(b_coefficients_[1]*sin(2 * PI*frequency) + b_coefficients_[2]*sin(4 * PI*frequency), 2)) /
		(pow((1 + a_coefficients_[1] * cos(2 * PI*frequency) + a_coefficients_[2] * cos(4 * PI*frequency)), 2) 
			+ pow(a_coefficients_[1] * sin(2 * PI*frequency) + a_coefficients_[2] * sin(4 * PI*frequency), 2)));
}

void IIRLowPassFilter::InitAcoefficients()
{
	if (a_coefficients_ == nullptr)
		a_coefficients_ = new double[COEFFICIENTS_NUMBER];

	
	
	a_coefficients_[0] = 1;
	a_coefficients_[1] = GetA1CoefficientValue();
	a_coefficients_[2] = GetA2Coefficient();
}

double IIRLowPassFilter::GetB0Coefficient() const
{
	return (angular_cutoff_frequency_ * angular_cutoff_frequency_) / GetCoefficientDenominator();
}

void IIRLowPassFilter::SetAbsorbtionFactor(double absorbiton_db)
{
	auto expression = pow(10, absorbiton_db / 10) - 1;
	absorbtion_factor_ = pow(expression, 1.0 / 4.0);
}

float IIRLowPassFilter::FilterOutput(float* memory, float sample) const
{
	
	auto tmp = (sample - a_coefficients_[1] * memory[0] - a_coefficients_[2] * memory[1]);
	auto ans= tmp*b_coefficients_[0] + b_coefficients_[1]*memory[0] + b_coefficients_[2]*memory[1];
	memory[1] = memory[0];
	memory[0] = float(tmp);

	return float(ans);
}

void IIRLowPassFilter::InitBcoefficients()
{
	if (b_coefficients_ == nullptr)
		b_coefficients_ = new double[MEMORY_SIZE];

	b_coefficients_[0] = GetB0Coefficient();
	b_coefficients_[1] = 2 * GetB0Coefficient();
	b_coefficients_[2] = GetB0Coefficient();
}

void IIRLowPassFilter::InitMemory()
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
