#include "digital_filter.h"


DigitalFilter::DigitalFilter(double cutoff_frequency, double sample_rate, bypassState bypass_state) : a_coefficients_(nullptr), 
																									  b_coefficients_(nullptr), 
																									  memory_left_(nullptr),
																									  memory_right_(nullptr),
																									  cutoff_frequency_(cutoff_frequency),
																									  sample_rate_(sample_rate),
																									  bypass_state_(bypass_state)			
{
	assert(sample_rate_ > 0);
	assert(cutoff_frequency_ > 0);
	SetAngularCutoffFrequency();
}

DigitalFilter::~DigitalFilter()
{
	if (a_coefficients_ != nullptr)
		delete a_coefficients_;
	if (b_coefficients_ != nullptr)
		delete b_coefficients_;
	if (memory_left_ != nullptr)
		delete memory_left_;
	if (memory_right_ != nullptr)
		delete memory_right_;
}

void DigitalFilter::SetBypassState(bypassState new_bypass_state)
{
	bypass_state_ = new_bypass_state;
}

bypassState DigitalFilter::GetBypassState() const
{
	return bypass_state_;
}

void DigitalFilter::SetAngularCutoffFrequency()
{
	assert(sample_rate_ > 0);
	assert(cutoff_frequency_ > 0);
	angular_cutoff_frequency_ = 2 * sample_rate_*std::tan(PI*cutoff_frequency_ / sample_rate_);
}

