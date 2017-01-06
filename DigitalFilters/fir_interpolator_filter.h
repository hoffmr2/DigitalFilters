#pragma once
#include "digital_filter.h"
class FIRInterpolatorFilter :
	public DigitalFilter
{
public:
	FIRInterpolatorFilter(int filter_size,int interpolation_factor);
	~FIRInterpolatorFilter();

	float FilterOutput(float * samples, int start, int samples_vector_size) const;

	void InitBcoefficients() override;
	void InitMemory() override;
	void ShiftMemory(float sample) const;

	void ChangeCutoffFrequency(double newFpass) override {};
	float FilterOutputLeft(float sample) override { return 0; };
	float FilterOutputRight(float sample) override { return 0; };
	double Spectrum(double frequency) override { return 0; };
	void InitAcoefficients() override {} ;

protected:

	int filter_size_;
	int interpolation_factor_;
};

