#include "iir_band_pass_filter.h"


IIRBandPassFilter::IIRBandPassFilter(double cutoff_frequency, double sample_rate, double gain, double q_factor, bypassState bypass_state)
	: IIRParametricBandPassFilter(cutoff_frequency,sample_rate,gain,q_factor,bypass_state)
{
	ChangeQFactor(q_factor);
}

float IIRBandPassFilter::FilterOutputLeft(float sample)
{
	if (bypass_state_ != off)
		return sample;
	double ans;
	auto tmp = (sample - a_coefficients_[1] * memory_left_[0] - a_coefficients_[2] * memory_left_[1]);
	ans = tmp*b_coefficients_[0] + b_coefficients_[1]*memory_left_[0] + b_coefficients_[2]*memory_left_[1];
	memory_left_[1] = memory_left_[0];
	memory_left_[0] = float(tmp);


	return float(ans);
}

float IIRBandPassFilter::FilterOutputRight(float sample)
{
	if (bypass_state_ != off)
		return sample;

	double ans;
	auto tmp = (sample - a_coefficients_[1] * memory_right_[0] - a_coefficients_[2] * memory_right_[1]);
	ans = tmp*b_coefficients_[0] + b_coefficients_[1] * memory_right_[0] + b_coefficients_[2] * memory_right_[1];
	memory_right_[1] = memory_right_[0];
	memory_right_[0] = float(tmp);


	return float(ans);
}

double IIRBandPassFilter::Spectrum(double frequency)
{
	if (bypass_state_ == off)
	{
		double b0 = b_coefficients_[0], b1 = b_coefficients_[1], b2 = b_coefficients_[2];

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

void IIRBandPassFilter::InitBcoefficients()
{
	if (b_coefficients_ == nullptr)
		b_coefficients_ = new double[COEFFICIENTS_NUMBER];

	auto denominator = 4 * sample_rate_*sample_rate_ + angular_cutoff_frequency_*angular_cutoff_frequency_ + 2 * bandwidth_*sample_rate_;
	b_coefficients_[0] = 2 * bandwidth_*(sample_rate_ / denominator) / denominator;
	b_coefficients_[1] = 0;
	b_coefficients_[2] = -2 * bandwidth_*(sample_rate_ / denominator) / denominator;
}


