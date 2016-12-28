#pragma once

#ifndef IIR_SHELF_FILTER_H_
#define IIR_SHELF_FILTER_H_



#include "digital_filter.h"
class IIRShelfFilter :
	public DigitalFilter
{
public:
	static const int COEFFICIENTS_NUMBER = 2, MEMORY_SIZE = 1;
	IIRShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state = off);
	virtual ~IIRShelfFilter();

	void SetGainDb(double gain_db);

	
	float FilterOutputLeft(float sample, double b0,double b1) const;
	float FilterOutputRight(float sample, double b0, double b1) const;
	void ChangeCutoffFrequency(double cutoff_frequency) override;
	double Spectrum(double frequency, double b0,double b1) const;
	void InitAcoefficients() override;
	void InitBcoefficients() override;
	void InitMemory() override;
protected:
	bool gain_flag_;
	double gain_;
};

#endif
