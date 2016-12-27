#include "iir_low_shelf_filter.h"



IIRLowShelfFilter::IIRLowShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state) : IIRShelfFilter(cutoff_frequency,sample_rate,gain_db,bypass_state)
{
}

IIRLowShelfFilter::~IIRLowShelfFilter()
{
}

/*
 * calculates b coefficients taking care of gain
 */
void IIRLowShelfFilter::CalculateCorrectedCoefficients(double& b0, double& b1) const
{
	b0 = b_coefficients_[0] + b_coefficients_[1]*gain_;
	b1 = b_coefficients_[1]*gain_ - b_coefficients_[0];
}


/*
 * Returns left sample filtered by low shelf filter
 */
float IIRLowShelfFilter::FilterOutputLeft(float sample)
{
	double b0,b1;

	CalculateCorrectedCoefficients(b0, b1);
	if (bypass_state_ == off)
		return IIRShelfFilter::FilterOutputLeft(sample, b0, b1);
	else
		return sample;
	
}


/*
* Returns right sample filtered by low shelf filter
*/
float IIRLowShelfFilter::FilterOutputRight(float sample)
{
	double b0,b1;
	CalculateCorrectedCoefficients(b0, b1);
	if (bypass_state_ == off)
		return IIRShelfFilter::FilterOutputRight(sample, b0, b1);
	else
		return sample;
}
/*
* Returns amplitude spectrum value
* for given frequency
*/
double IIRLowShelfFilter::Spectrum(double frequency)
{
	double b0,b1;
	CalculateCorrectedCoefficients(b0, b1);
	if (bypass_state_ == off)
		return IIRShelfFilter::Spectrum(frequency, b0, b1);
	else
		return 1;
}
