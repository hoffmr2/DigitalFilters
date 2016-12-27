#include "iir_shelf_filter.h"


IIRShelfFilter::IIRShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state) : DigitalFilter(cutoff_frequency,sample_rate,bypass_state), gain_flag_(true)																														
{
	SetGainDb(gain_db);
	InitAcoefficients();
	InitBcoefficients();
	InitMemory();

}

IIRShelfFilter::~IIRShelfFilter()
{
}

void IIRShelfFilter::SetGainDb(double gain_db)
{
	
	gain_flag_ = (gain_db <= 0) ? true : false;
	auto exp = (gain_db) / 20.0;
	gain_ = std::pow(10, (-std::abs(exp)));
}

float IIRShelfFilter::FilterOutputLeft(float sample, double b0, double b1) const
{
	double ans;
	if (gain_flag_ == true)
	{
		ans = (sample - a_coefficients_[1]*memory_left_[0])*b0 + b1*memory_left_[0];
		memory_left_[0] = sample - a_coefficients_[1]*memory_left_[0];
	}
	else
	{
		ans = (sample - b1*memory_left_[0] / b0) / b0 + a_coefficients_[1]*memory_left_[0] / b0;
		memory_left_[0] = sample - b1*memory_left_[0] / b0;
	}
	return ans;
}

void IIRShelfFilter::ChangeCutoffFrequency(double cutoff_frequency)
{
	assert(cutoff_frequency > 0);
	cutoff_frequency_ = cutoff_frequency;
	SetAngularCutoffFrequency();
	InitAcoefficients();
	InitBcoefficients();
}

void IIRShelfFilter::InitAcoefficients()
{
	if (a_coefficients_ == nullptr)
		a_coefficients_ = new double[COEFFICIENTS_NUMBER];
	
	a_coefficients_[0] = 1;
	a_coefficients_[1] =  (angular_cutoff_frequency_ - 2 * sample_rate_) / (angular_cutoff_frequency_ + 2 * sample_rate_);
}

void IIRShelfFilter::InitBcoefficients()
{
	if(b_coefficients_ == nullptr)
		b_coefficients_ = new double[COEFFICIENTS_NUMBER];
	b_coefficients_[0] =  2 * sample_rate_ / (angular_cutoff_frequency_ + 2 * sample_rate_);
	b_coefficients_[1] = 1- b_coefficients_[0];

}

void IIRShelfFilter::InitMemory()
{
	if(memory_left_ == nullptr)
		memory_left_ = new float[MEMORY_SIZE];
	if(memory_left_ ==nullptr)
		memory_right_ = new float[MEMORY_SIZE];
	memory_left_[0] = 0;
	memory_right_[0] = 0;
}
