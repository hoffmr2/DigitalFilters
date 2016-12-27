#pragma once

#ifndef IIR_HIGH_SHELF_FILTER_H_
#define IIR_HIGH_SHELF_FILTER_H_


#include "iir_shelf_filter.h"
class IIRHighShelfFilter :
	public IIRShelfFilter
{
public:
	IIRHighShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state = off);
	virtual ~IIRHighShelfFilter();
	void CalculateCorrectedCoefficients(double& b0, double& b1) const;
	float FilterOutputLeft(float sample) override;
	float FilterOutputRight(float sample) override;
	double Spectrum(double frequency) override;
};

#endif
