#include "fir_low_pass_filter.h"


void FirLowPassFilter::InitFilter()
{
	InitDFactor();
	CalculateFilterSize();
	InitMemory();
	InitBetaFactor();
	InitBcoefficients();
}

FirLowPassFilter::FirLowPassFilter(double sample_rate, double pass_band_frequency, double stop_band_frequency, double absorbtion_in_stop_band)
	: FirFilter(sample_rate, pass_band_frequency), stop_band_frequency_(stop_band_frequency), absorbtion_in_stop_band_(absorbtion_in_stop_band)
{
	InitFilter();

}

FirLowPassFilter::~FirLowPassFilter()
{
}

void FirLowPassFilter::InitBcoefficients()
{
	double h, w;

	if (b_coefficients_ != nullptr)
		delete b_coefficients_;
	b_coefficients_ = new double[filter_size_];
	auto fc = (cutoff_frequency_ / 2 + stop_band_frequency_ / 2) / sample_rate_;
	double half_of_filter_size = filter_size_ / 2;
	auto besseli_value = BesselZeroKindFunction(beta_factor_);
	for (int i = 0; i < filter_size_; ++i)
	{
		if (i - half_of_filter_size == 0)
			h = 2 * fc;
		else
			h = 2 * fc*sin(2 * 3.14*fc*(i - half_of_filter_size)) / (2 * 3.14*fc*(i - half_of_filter_size));
		w = BesselZeroKindFunction(beta_factor_*sqrt(1 - pow((i - half_of_filter_size) / half_of_filter_size, 2))) / besseli_value;
		b_coefficients_[i] = (w*h);
	}
}

void FirLowPassFilter::ChangeCutoffFrequency(double passband_requency, double stopband_frequency)
{
	cutoff_frequency_ = passband_requency;
	stop_band_frequency_ = stopband_frequency;
	InitFilter();
}

double FirLowPassFilter::factorial(unsigned arg)
{
	if (arg == 0 || arg == 1) return 1;
	int ans = 2;
	for (unsigned int i = 3; i <= arg; ++i)
		ans *= i;
	return ans;
}

void FirLowPassFilter::ChangeCutoffFrequency(double newFpass)
{

}

void FirLowPassFilter::CalculateFilterSize()
{
	assert(d_factor_ != 0);
	assert(stop_band_frequency_ > cutoff_frequency_);
	auto df = stop_band_frequency_ - cutoff_frequency_;
	filter_size_ = static_cast<int>(ceil((d_factor_ * sample_rate_) / (df)));
	filter_size_ += (filter_size_ % 2 == 0) ? 1 : 0;
}

void FirLowPassFilter::InitDFactor()
{
	assert(absorbtion_in_stop_band_ != 0);
	d_factor_ = (absorbtion_in_stop_band_ > 21) ? (absorbtion_in_stop_band_ - 7.95) / 14.36 : 0.922;
}

void FirLowPassFilter::InitBetaFactor()
{
	assert(absorbtion_in_stop_band_ != 0);
	if (absorbtion_in_stop_band_ < 21)
	{
		beta_factor_ = 0.0;
		return;
	}
	if (absorbtion_in_stop_band_ >= 21 && absorbtion_in_stop_band_ <= 51)
	{
		beta_factor_ = 0.5842*pow(absorbtion_in_stop_band_ - 21, 0.4) + 0.07886*(absorbtion_in_stop_band_ - 21);
		return;
	}
	beta_factor_ = 0.1102*(absorbtion_in_stop_band_ - 8.7);
}

double FirLowPassFilter::BesselZeroKindFunction(double beta) const
{
	double ans = 1;
	for (int k = 1; k < 7; ++k)
	{
		ans += pow(pow(beta / 2, k) / factorial(k), 2);
	}
	return ans;
}