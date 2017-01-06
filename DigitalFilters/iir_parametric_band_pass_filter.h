#pragma once

#ifndef IIR_PARAMETRIC_BAND_PASS_FILTER_H_
#define IIR_PARAMETRIC_BAND_PASS_FILTER_H_



#include "digital_filter.h"
class IIRParametricBandPassFilter :
	public DigitalFilter
{
public:
	static const int COEFFICIENTS_NUMBER = 3, MEMORY_SIZE= 2;
	IIRParametricBandPassFilter(double cutoff_frequency, double sample_rate, double gain, double q_factor, bypassState bypass_state = off);
	virtual ~IIRParametricBandPassFilter();
	void CalculateBandwidth();

	void ChangeCutoffFrequency(double cutoff_frequency) override;
	void ChangeQFactor(double q_factor);
	float FilterOutputLeft(float sample) override;
	void CalculateCorrectedBCoefficients(double& b0, double& b1, double& b2) const;
	float FilterOutputRight(float sample) override;
	float FilterOutput(float sample, float* memory) const;
	double Spectrum(double frequency) override;
	void InitAcoefficients() override;
	void InitBcoefficients() override;
	void InitMemory() override;
	void SetGainDb(double gain_db);

protected:
	double q_factor_;
	double bandwidth_;
	bool gain_flag_;
	double gain_;
};

#endif
