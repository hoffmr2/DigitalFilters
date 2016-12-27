#pragma once

/*
 * Digital Filter is a base abstract base class
 * for creating digital filters
 */

#ifndef DIGITAL_FILTER_H_
#define DIGITAL_FILTER_H_

#define PI 3.141592653589793
enum bypassState { on, off };


#include <assert.h>
#include <cmath>

class DigitalFilter
{
public:
	DigitalFilter(double cutoff_frequency, double sample_rate, bypassState bypass_state = off);
	virtual ~DigitalFilter();

	void SetBypassState(bypassState new_bypass_state);
	bypassState GetBypassState() const;



	virtual void ChangeCutoffFrequency(double newFpass)=0;
	virtual float FilterOutputLeft(float sample) = 0;
	virtual float FilterOutputRight(float sample) = 0;
	virtual double Spectrum(double frequency) = 0;
	virtual void InitAcoefficients() = 0;
	virtual void InitBcoefficients() = 0;
	virtual void InitMemory() = 0;

protected:
	virtual void SetAngularCutoffFrequency();

	//Private variables
	double* a_coefficients_;
	double* b_coefficients_;
	float* memory_left_;
	float* memory_right_;

	double cutoff_frequency_;
	double sample_rate_;
	double angular_cutoff_frequency_;

	bypassState bypass_state_;

};

#endif
