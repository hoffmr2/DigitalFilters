#pragma once
#ifndef FIR_BAND_PASS_FILTER_H_
#define FIR_BAND_PASS_FILTER_H_

#include "fir_filter.h"
class FirBandPassFilter :
  public FirFilter
{
public:
  void InitFilter();
  FirBandPassFilter(double sample_rate, double pass_band_frequency, double stop_band_frequency, double absorbtion_in_stop_band, double middle_frequency);
  ~FirBandPassFilter();

  virtual void InitBcoefficients() override;
  void ChangeCutoffFrequency(double passband_requency, double stopband_frequency,double middle_frequency);

  void CalculateFilterSize();
  void InitDFactor();
  void InitBetaFactor();
protected:
  double stop_band_frequency_;
  double absorbtion_in_stop_band_;
  double middle_frequency_;
  double beta_factor_;
  double d_factor_;
private:
  static double factorial(unsigned int arg);
  void ChangeCutoffFrequency(double newFpass) override;
  double BesselZeroKindFunction(double beta) const;
};

#endif