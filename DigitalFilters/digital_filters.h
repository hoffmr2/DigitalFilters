#pragma once
#ifndef DIGITAL_FILTERS_H_
#define DIGITALFILTERS_H_

#define PI 3.141592653589793
#include <cassert>

namespace HoffFilters
{
	enum bypassState { on, off };


#include <assert.h>
#include <cmath>
	using namespace std;

	class DigitalFilter
	{
	public:
		DigitalFilter(double cutoff_frequency, double sample_rate, bypassState bypass_state = off);
		virtual ~DigitalFilter();

		void SetBypassState(bypassState new_bypass_state);
		bypassState GetBypassState() const;



		virtual void ChangeCutoffFrequency(double newFpass) = 0;
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


	DigitalFilter::DigitalFilter(double cutoff_frequency, double sample_rate, bypassState bypass_state) : a_coefficients_(nullptr),
		b_coefficients_(nullptr),
		memory_left_(nullptr),
		memory_right_(nullptr),
		cutoff_frequency_(cutoff_frequency),
		sample_rate_(sample_rate),
		bypass_state_(bypass_state)
	{
		assert(sample_rate_ > 0);
		assert(cutoff_frequency_ > 0);
		SetAngularCutoffFrequency();
	}

	DigitalFilter::~DigitalFilter()
	{
		if (a_coefficients_ != nullptr)
			delete a_coefficients_;
		if (b_coefficients_ != nullptr)
			delete b_coefficients_;
		if (memory_left_ != nullptr)
			delete memory_left_;
		if (memory_right_ != nullptr)
			delete memory_right_;
	}

	void DigitalFilter::SetBypassState(bypassState new_bypass_state)
	{
		bypass_state_ = new_bypass_state;
	}

	bypassState DigitalFilter::GetBypassState() const
	{
		return bypass_state_;
	}

	void DigitalFilter::SetAngularCutoffFrequency()
	{
		assert(sample_rate_ > 0);
		assert(cutoff_frequency_ > 0);
		angular_cutoff_frequency_ = 2 * sample_rate_*tan(PI*cutoff_frequency_ / sample_rate_);
	}



	class IIRShelfFilter :
		public DigitalFilter
	{
	public:
		static const int COEFFICIENTS_NUMBER = 2, MEMORY_SIZE = 1;
		IIRShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state = off);
		virtual ~IIRShelfFilter();

		void SetGainDb(double gain_db);


		float FilterOutputLeft(float sample, double b0, double b1) const;
		float FilterOutputRight(float sample, double b0, double b1) const;
		void ChangeCutoffFrequency(double cutoff_frequency) override;
		double Spectrum(double frequency, double b0, double b1) const;
		void InitAcoefficients() override;
		void InitBcoefficients() override;
		void InitMemory() override;
	protected:
		bool gain_flag_;
		double gain_;
	};


	IIRShelfFilter::IIRShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state) : DigitalFilter(cutoff_frequency, sample_rate, bypass_state), gain_flag_(true)
	{
		SetGainDb(gain_db);
		InitAcoefficients();
		InitBcoefficients();
		InitMemory();

	}

	IIRShelfFilter::~IIRShelfFilter()
	{
	}

	/*
	* Sets filter gain
	* arguments must be given in [DB]
	*/
	void IIRShelfFilter::SetGainDb(double gain_db)
	{

		gain_flag_ = (gain_db <= 0) ? true : false;
		auto exp = (gain_db) / 20.0;
		gain_ = pow(10, (-abs(exp)));
	}

	/*
	* Returns the value od filtered sample
	* for left channel
	*/
	float IIRShelfFilter::FilterOutputLeft(float sample, double b0, double b1) const
	{
		double ans;
		if (gain_flag_ == true)
		{
			ans = (sample - a_coefficients_[1] * memory_left_[0])*b0 + b1*memory_left_[0];
			memory_left_[0] = float(sample - a_coefficients_[1] * memory_left_[0]);
		}
		else
		{
			ans = (sample - b1*memory_left_[0] / b0) / b0 + a_coefficients_[1] * memory_left_[0] / b0;
			memory_left_[0] = float(sample - b1*memory_left_[0] / b0);
		}
		return float(ans);
	}

	/*
	* returns the value of filtered sample
	* for right channel
	*/
	float IIRShelfFilter::FilterOutputRight(float sample, double b0, double b1) const
	{
		double ans;
		if (gain_flag_ == true)
		{
			ans = (sample - a_coefficients_[1] * memory_right_[0])*b0 + b1*memory_right_[0];
			memory_right_[0] = float(sample - a_coefficients_[1] * memory_right_[0]);
		}
		else
		{
			ans = (sample - b1*memory_right_[0] / b0) / b0 + a_coefficients_[1] * memory_right_[0] / b0;
			memory_right_[0] = float(sample - b1*memory_right_[0] / b0);
		}
		return float(ans);
	}


	/*
	* Changes cutoff frequency to new one
	* and sets all necessary parameters
	*/
	void IIRShelfFilter::ChangeCutoffFrequency(double cutoff_frequency)
	{
		assert(cutoff_frequency > 0);
		cutoff_frequency_ = cutoff_frequency;
		SetAngularCutoffFrequency();
		InitAcoefficients();
		InitBcoefficients();
	}

	/*
	* returns the amplitude spectrum value
	* for given frequency and b0 and b1 filter coefficients
	*/
	double IIRShelfFilter::Spectrum(double frequency, double b0, double b1) const
	{
		if (gain_flag_ == true)
			return sqrt((pow(b0 + b1 *cos(2 * PI*frequency), 2) + pow(-b1 *sin(2 * PI*frequency), 2))
				/ (pow(a_coefficients_[0] + a_coefficients_[1] * cos(2 * PI*frequency), 2) + pow(-a_coefficients_[1] * sin(2 * PI*frequency), 2)));
		else
			return sqrt((pow(a_coefficients_[0] + a_coefficients_[1] * cos(2 * PI*frequency), 2) + pow(-a_coefficients_[1] * sin(2 * PI*frequency), 2))
				/ (pow(b0 + b1 *cos(2 * PI*frequency), 2) + pow(-b1 *sin(2 * PI*frequency), 2)));
	}

	void IIRShelfFilter::InitAcoefficients()
	{
		if (a_coefficients_ == nullptr)
			a_coefficients_ = new double[COEFFICIENTS_NUMBER];

		a_coefficients_[0] = 1;
		a_coefficients_[1] = (angular_cutoff_frequency_ - 2 * sample_rate_) / (angular_cutoff_frequency_ + 2 * sample_rate_);
	}

	void IIRShelfFilter::InitBcoefficients()
	{
		if (b_coefficients_ == nullptr)
			b_coefficients_ = new double[COEFFICIENTS_NUMBER];
		b_coefficients_[0] = 2 * sample_rate_ / (angular_cutoff_frequency_ + 2 * sample_rate_);
		b_coefficients_[1] = 1 - b_coefficients_[0];

	}

	void IIRShelfFilter::InitMemory()
	{
		if (memory_left_ == nullptr)
			memory_left_ = new float[MEMORY_SIZE];
		if (memory_right_ == nullptr)
			memory_right_ = new float[MEMORY_SIZE];
		memory_left_[0] = 0;
		memory_right_[0] = 0;
	}


	class IIRLowShelfFilter :
		public IIRShelfFilter
	{
	public:
		IIRLowShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state = off);
		virtual ~IIRLowShelfFilter();
		void CalculateCorrectedCoefficients(double& b0, double& b1) const;
		float FilterOutputLeft(float sample) override;
		float FilterOutputRight(float sample) override;
		double Spectrum(double frequency) override;
	};

	IIRLowShelfFilter::IIRLowShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state) : IIRShelfFilter(cutoff_frequency, sample_rate, gain_db, bypass_state)
	{
	}

	IIRLowShelfFilter::~IIRLowShelfFilter()
	{
	}

	/*
	* calculates b coefficients taking care of gain
	*/
	void IIRLowShelfFilter::CalculateCorrectedCoefficients(double& b0, double& b1) const
	{
		b0 = b_coefficients_[0] + b_coefficients_[1] * gain_;
		b1 = b_coefficients_[1] * gain_ - b_coefficients_[0];
	}


	/*
	* Returns left sample filtered by low shelf filter
	*/
	float IIRLowShelfFilter::FilterOutputLeft(float sample)
	{
		double b0, b1;

		CalculateCorrectedCoefficients(b0, b1);
		if (bypass_state_ == off)
			return IIRShelfFilter::FilterOutputLeft(sample, b0, b1);
		else
			return sample;

	}


	/*
	* Returns right sample filtered by low shelf filter
	*/
	float IIRLowShelfFilter::FilterOutputRight(float sample)
	{
		double b0, b1;
		CalculateCorrectedCoefficients(b0, b1);
		if (bypass_state_ == off)
			return IIRShelfFilter::FilterOutputRight(sample, b0, b1);
		else
			return sample;
	}
	/*
	* Returns amplitude spectrum value
	* for given frequency
	*/
	double IIRLowShelfFilter::Spectrum(double frequency)
	{
		double b0, b1;
		CalculateCorrectedCoefficients(b0, b1);
		if (bypass_state_ == off)
			return IIRShelfFilter::Spectrum(frequency, b0, b1);
		else
			return 1;
	}


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


	IIRHighShelfFilter::IIRHighShelfFilter(double cutoff_frequency, double sample_rate, double gain_db, bypassState bypass_state) : IIRShelfFilter(cutoff_frequency, sample_rate, gain_db, bypass_state)
	{
	}

	IIRHighShelfFilter::~IIRHighShelfFilter()
	{
	}

	/*
	* calculates b coefficients taking care of gain
	*/
	void IIRHighShelfFilter::CalculateCorrectedCoefficients(double& b0, double& b1) const
	{
		b0 = b_coefficients_[0] * gain_ + b_coefficients_[1];
		b1 = b_coefficients_[1] - b_coefficients_[0] * gain_;
	}


	/*
	* Returns left sample filtered by low shelf filter
	*/
	float IIRHighShelfFilter::FilterOutputLeft(float sample)
	{
		double b0, b1;

		CalculateCorrectedCoefficients(b0, b1);
		if (bypass_state_ == off)
			return IIRShelfFilter::FilterOutputLeft(sample, b0, b1);
		else
			return sample;

	}


	/*
	* Returns right sample filtered by low shelf filter
	*/
	float IIRHighShelfFilter::FilterOutputRight(float sample)
	{
		double b0, b1;
		CalculateCorrectedCoefficients(b0, b1);
		if (bypass_state_ == off)
			return IIRShelfFilter::FilterOutputRight(sample, b0, b1);
		else
			return sample;
	}
	/*
	* Returns amplitude spectrum value
	* for given frequency
	*/
	double IIRHighShelfFilter::Spectrum(double frequency)
	{
		double b0, b1;
		CalculateCorrectedCoefficients(b0, b1);
		if (bypass_state_ == off)
			return IIRShelfFilter::Spectrum(frequency, b0, b1);
		else
			return 1;
	}


	class IIRParametricBandPassFilter :
		public DigitalFilter
	{
	public:
		static const int COEFFICIENTS_NUMBER = 3, MEMORY_SIZE = 2;
		IIRParametricBandPassFilter(double cutoff_frequency, double sample_rate, double gain, double q_factor, bypassState bypass_state = off);
		virtual ~IIRParametricBandPassFilter();
		void CalculateBandwidth();

		void ChangeCutoffFrequency(double cutoff_frequency) override;
		void ChangeQFactor(double q_factor);
		float FilterOutputLeft(float sample) override;
		void CalculateCorrectedBCoefficients(double& b0, double& b1, double& b2) const;
		float FilterOutputRight(float sample) override;
		float FilterOutput(float sample, float* memory);
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


	IIRParametricBandPassFilter::IIRParametricBandPassFilter(double cutoff_frequency, double sample_rate, double gain, double q_factor, bypassState bypass_state) : DigitalFilter(cutoff_frequency, sample_rate, bypass_state),
		q_factor_(q_factor),
		gain_flag_(true)
	{
		SetAngularCutoffFrequency();
		CalculateBandwidth();
		SetGainDb(gain);
		InitMemory();
		InitAcoefficients();
		InitBcoefficients();
	}

	IIRParametricBandPassFilter::~IIRParametricBandPassFilter()
	{
	}

	void IIRParametricBandPassFilter::CalculateBandwidth()
	{
		assert(q_factor_ >= 0);
		bandwidth_ = q_factor_*angular_cutoff_frequency_;
	}

	void IIRParametricBandPassFilter::ChangeCutoffFrequency(double cutoff_frequency)
	{
		assert(cutoff_frequency > 0);
		cutoff_frequency_ = cutoff_frequency;

		SetAngularCutoffFrequency();
		CalculateBandwidth();
		InitAcoefficients();
		InitBcoefficients();
	}

	void IIRParametricBandPassFilter::ChangeQFactor(double q_factor)
	{
		assert(q_factor >= 0);
		q_factor_ = q_factor;

		SetAngularCutoffFrequency();
		CalculateBandwidth();
		InitAcoefficients();
		InitBcoefficients();
	}

	float IIRParametricBandPassFilter::FilterOutputLeft(float sample)
	{
		if (bypass_state_ == off)
			return FilterOutput(sample, memory_left_);
		else
			return sample;

	}

	void IIRParametricBandPassFilter::CalculateCorrectedBCoefficients(double& b0, double& b1, double& b2) const
	{
		b0 = b_coefficients_[0] + gain_*b_coefficients_[1];
		b1 = b_coefficients_[2];
		b2 = b_coefficients_[0] - gain_*b_coefficients_[1];
	}

	float IIRParametricBandPassFilter::FilterOutputRight(float sample)
	{
		if (bypass_state_ == off)
			return FilterOutput(sample, memory_right_);
		else
			return sample;
	}

	float IIRParametricBandPassFilter::FilterOutput(float sample, float* memory)
	{
		double ans;
		double b0, b1, b2;
		CalculateCorrectedBCoefficients(b0, b1, b2);


		if (gain_flag_ == true)
		{
			auto tmp = (sample - a_coefficients_[1] * memory[0] - a_coefficients_[2] * memory[1]);
			ans = tmp*b0 + b1*memory[0] + b2*memory[1];
			memory[1] = memory[0];
			memory[0] = float(tmp);
		}
		else
		{
			auto tmp = (sample - b1 * memory[0] / b0 - b2 * memory[1] / b0);
			ans = tmp / b0 + a_coefficients_[1] * memory[0] / b0 + a_coefficients_[2] * memory[1] / b0;
			memory[1] = memory[0];
			memory[0] = float(tmp);
		}

		return float(ans);
	}

	double IIRParametricBandPassFilter::Spectrum(double frequency)
	{
		if (bypass_state_ == off)
		{
			double b0, b1, b2;
			CalculateCorrectedBCoefficients(b0, b1, b2);

			if (gain_flag_ == true)
				return sqrt((pow((b0 + b1*cos(2 * PI*frequency) + b2*cos(4 * PI*frequency)), 2) + pow(b1*sin(2 * PI*frequency) + b2*sin(4 * PI*frequency), 2)) /
				(pow((1 + a_coefficients_[1] * cos(2 * PI*frequency) + a_coefficients_[2] * cos(4 * PI*frequency)), 2) + pow(a_coefficients_[1] * sin(2 * PI*frequency) + a_coefficients_[2] * sin(4 * PI*frequency), 2)));
			else
				return sqrt((pow((1 + a_coefficients_[1] * cos(2 * PI*frequency) + a_coefficients_[2] * cos(4 * PI*frequency)), 2) + pow(a_coefficients_[1] * sin(2 * PI*frequency) + a_coefficients_[2] * sin(4 * PI*frequency), 2)) /
				(pow((b0 + b1*cos(2 * PI*frequency) + b2*cos(4 * PI*frequency)), 2) + pow(b1*sin(2 * PI*frequency) + b2*sin(4 * PI*frequency), 2)));
		}
		else
			return 1;
	}

	void IIRParametricBandPassFilter::InitAcoefficients()
	{

		if (a_coefficients_ == nullptr)
			a_coefficients_ = new double[COEFFICIENTS_NUMBER];

		auto denominator = 4 * sample_rate_*sample_rate_ + angular_cutoff_frequency_*angular_cutoff_frequency_ + sample_rate_*(2 * bandwidth_);
		a_coefficients_[0] = 1;
		a_coefficients_[1] = (angular_cutoff_frequency_*(2 * angular_cutoff_frequency_) / denominator) - sample_rate_*((8 * sample_rate_) / denominator);
		a_coefficients_[2] = ((4 * sample_rate_) / denominator)*sample_rate_ + angular_cutoff_frequency_*(angular_cutoff_frequency_ / denominator) - bandwidth_*((2 * sample_rate_) / denominator);
	}

	void IIRParametricBandPassFilter::InitBcoefficients()
	{
		if (b_coefficients_ == nullptr)
			b_coefficients_ = new double[COEFFICIENTS_NUMBER];

		auto denominator = 4 * sample_rate_*sample_rate_ + angular_cutoff_frequency_*angular_cutoff_frequency_ + 2 * bandwidth_*sample_rate_;
		b_coefficients_[0] = ((4 * sample_rate_*sample_rate_) / denominator) + ((angular_cutoff_frequency_*angular_cutoff_frequency_) / denominator);
		b_coefficients_[1] = 2 * bandwidth_*(sample_rate_ / denominator);
		b_coefficients_[2] = ((2 * angular_cutoff_frequency_*angular_cutoff_frequency_) / denominator) - ((8 * sample_rate_*sample_rate_) / denominator);
	}

	void IIRParametricBandPassFilter::InitMemory()
	{
		if (memory_left_ == nullptr)
			memory_left_ = new float[MEMORY_SIZE];
		if (memory_right_ == nullptr)
			memory_right_ = new float[MEMORY_SIZE];
		for (auto i = 0; i < MEMORY_SIZE; ++i)
		{
			memory_left_[i] = 0;
			memory_right_[i] = 0;
		}
	}

	void IIRParametricBandPassFilter::SetGainDb(double gain_db)
	{
		gain_flag_ = (gain_db <= 0) ? true : false;
		auto exp = (gain_db) / 20.0;
		gain_ = pow(10, (-abs(exp)));
	}




}

#endif
