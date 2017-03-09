#pragma once
#ifndef IIR_BAND_PASS_FILTER_H_
#define IIR_BAND_PASS_FILTER_H_

#include "iir_parametric_band_pass_filter.h"
class IIRBandPassFilter :
	public IIRParametricBandPassFilter
{
public:
	IIRBandPassFilter(double cutoff_frequency, double sample_rate, double gain, double q_factor, bypassState bypass_state = off);

	float FilterOutputLeft(float sample) override;
	float FilterOutputRight(float sample) override;
	double Spectrum(double frequency) override;
	void InitBcoefficients() override;
};

#endif
