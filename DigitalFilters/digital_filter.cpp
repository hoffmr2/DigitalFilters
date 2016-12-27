#include "digital_filter.h"


DigitalFilter::DigitalFilter(double cutoff_frequency, double sample_rate, bypassState bypass_state) : cutoff_frequency_(cutoff_frequency), 
																									  sample_rate_(sample_rate), 
																									  bypass_state_(bypass_state)
{
	assert(sample_rate_ > 0);
	assert(cutoff_frequency_ > 0);
	SetAngularCutoffFrequency();
}

DigitalFilter::~DigitalFilter()
{
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

