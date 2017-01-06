#pragma once

#ifndef IIR_LOW_PASS_FILTER_H_
#define IIR_LOW_PASS_FILTER_H_

#include "digital_filter.h"
class IIRLowPassFilter :
	public DigitalFilter
{
public:
	
	IIRLowPassFilter(double cutoff_frequency, double sample_rate, double absorbtion_db = 3, bypassState bypass_state = off);
	~IIRLowPassFilter();

	virtual void ChangeCutoffFrequency(double newFpass) override;
	virtual float FilterOutputLeft(float sample) override;
	virtual float FilterOutputRight(float sample) override;
	virtual double Spectrum(double frequency) override;
	virtual void InitAcoefficients() override;
	
	virtual void InitBcoefficients() override;
	virtual void InitMemory() override;
protected:
	double GetCoefficientDenominator() const;
	double GetA1CoefficientValue() const;
	int GetA2Coefficient() const;
	double GetB0Coefficient() const;
	void SetAbsorbtionFactor(double absorbtion_db);
	float FilterOutput(float* memory,float sample) const;

	static const int MEMORY_SIZE = 2, COEFFICIENTS_NUMBER = 3;
	double absorbtion_factor_;
};

#endif
