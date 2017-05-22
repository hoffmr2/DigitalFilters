#include "fir_band_pass_filter.h"



void FirBandPassFilter::InitFilter()
{
  InitDFactor();
  CalculateFilterSize();
  InitMemory();
  InitBetaFactor();
  InitBcoefficients();
}

FirBandPassFilter::FirBandPassFilter(double sample_rate, double pass_band_frequency, double stop_band_frequency, double absorbtion_in_stop_band,double middle_frequency)
  : FirFilter(sample_rate, pass_band_frequency), stop_band_frequency_(stop_band_frequency), absorbtion_in_stop_band_(absorbtion_in_stop_band),middle_frequency_(middle_frequency)
{
  InitFilter();

}

FirBandPassFilter::~FirBandPassFilter()
{
  
}


void FirBandPassFilter::InitBcoefficients()
{
  double h, w,f0;

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
    f0 = std::cos(2 * PI*middle_frequency_*i/sample_rate_);
    b_coefficients_[i] = (w*h*f0);
  }
}

void FirBandPassFilter::ChangeCutoffFrequency(double passband_requency, double stopband_frequency,double middle_frequency)
{
  cutoff_frequency_ = passband_requency;
  stop_band_frequency_ = stopband_frequency;
  middle_frequency_ = middle_frequency;
  InitFilter();
}

double FirBandPassFilter::factorial(unsigned arg)
{
  if (arg == 0 || arg == 1) return 1;
  int ans = 2;
  for (unsigned int i = 3; i <= arg; ++i)
    ans *= i;
  return ans;
}

void FirBandPassFilter::ChangeCutoffFrequency(double newFpass)
{
  cutoff_frequency_ = newFpass;
  stop_band_frequency_ = 2 * cutoff_frequency_;
  InitFilter();
}

void FirBandPassFilter::CalculateFilterSize()
{
  assert(d_factor_ != 0);
  assert(stop_band_frequency_ > cutoff_frequency_);
  auto df = stop_band_frequency_ - cutoff_frequency_;
  filter_size_ = static_cast<int>(ceil((d_factor_ * sample_rate_) / (df)));
  filter_size_ += (filter_size_ % 2 == 0) ? 1 : 0;
}

void FirBandPassFilter::InitDFactor()
{
  assert(absorbtion_in_stop_band_ != 0);
  d_factor_ = (absorbtion_in_stop_band_ > 21) ? (absorbtion_in_stop_band_ - 7.95) / 14.36 : 0.922;
}

void FirBandPassFilter::InitBetaFactor()
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

double FirBandPassFilter::BesselZeroKindFunction(double beta) const
{
  double ans = 1;
  for (int k = 1; k < 7; ++k)
  {
    ans += pow(pow(beta / 2, k) / factorial(k), 2);
  }
  return ans;
}
