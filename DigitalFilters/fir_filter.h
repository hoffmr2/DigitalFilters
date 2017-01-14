#pragma once

#ifndef FIR_FILTER_H_
#define FIR_FILTER_H_

#define A_COEFFICIENTS_VALUE 1;


#include "digital_filter.h"
class FirFilter :
	public DigitalFilter
{
public:
	FirFilter(double sample_rate, double pass_band_frequency);
	FirFilter(double sample_rate, double pass_band_frequency, double* b_coefficients,int filter_size);
	~FirFilter();

	virtual float FilterOutputLeft(float sample) override;
	virtual float FilterOutputRight(float sample) override;
	virtual double Spectrum(double frequency) override;
	
	void InitBcoefficients() override;
	void InitAcoefficients() override;

	void ShiftMemory(float* memory) const;
	void ClearMemory() const;
	virtual void InitMemory() override;
protected:
	float FilterOutput(float sample, float* memory) const;

	int filter_size_;

private:
	void ChangeCutoffFrequency(double newFpass) override;
	double* tmp_b_coeffcients_;








};

#endif
