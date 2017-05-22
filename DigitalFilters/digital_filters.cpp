#include "digital_filters.h"

namespace HoffFilters
{
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

	FIRInterpolatorFilter::FIRInterpolatorFilter(int filter_size, int interpolation_factor) : DigitalFilter(1, 1), filter_size_(filter_size), interpolation_factor_(interpolation_factor)
	{
		InitMemory();
		InitBcoefficients();
	}

	FIRInterpolatorFilter::~FIRInterpolatorFilter()
	{
		if (b_coefficients_ != nullptr)
			delete[] b_coefficients_;
		if (memory_left_ != nullptr)
			delete[] memory_left_;
	}

	float FIRInterpolatorFilter::FilterOutput(float* samples, int start, int samples_vector_size) const
	{
		float ans = 0;
		for (int i = 0; i < filter_size_; ++i)
			ans += memory_left_[i] * b_coefficients_[i];
		ans += samples[start];
		ShiftMemory(samples[start]);
		for (int i = start + 1; i < start + filter_size_ && i < samples_vector_size; ++i)
			ans += samples[i] * b_coefficients_[i - start];

		return ans;
	}

	void FIRInterpolatorFilter::InitBcoefficients()
	{
		b_coefficients_ = new double[2 * filter_size_ + 1];

		for (auto i = 0; i < 2 * filter_size_ + 1; ++i)
		{
			if (i - filter_size_ != 0)
				b_coefficients_[i] = sin(PI*(i - filter_size_) / interpolation_factor_) / (PI*(i - filter_size_) / interpolation_factor_);
			else
				b_coefficients_[i] = 1;
			auto hanning_window = (0.54 - 0.46* cos(2 * PI*(i - filter_size_) / (2 * filter_size_)));
			b_coefficients_[i] *= hanning_window;
		}
	}

	void FIRInterpolatorFilter::InitMemory()
	{
		if (memory_left_ == nullptr)
			memory_left_ = new float[filter_size_];
		for (auto i = 0; i < filter_size_; ++i)
			memory_left_[i] = 0;
	}

	void FIRInterpolatorFilter::ShiftMemory(float sample) const
	{
		float tmp;
		for (auto i = filter_size_ - 2; i >= 0; --i)
		{
			tmp = memory_left_[i];
			memory_left_[i] = memory_left_[i + 1];
		}
		memory_left_[filter_size_ - 1] = sample;
	}

	IIRLowPassFilter::IIRLowPassFilter(double cutoff_frequency, double sample_rate, double absorbtion_db, bypassState bypass_state)
		:DigitalFilter(cutoff_frequency, sample_rate, bypass_state)
	{
		assert(absorbtion_db > 0);
		SetAbsorbtionFactor(absorbtion_db);
		InitMemory();
		InitAcoefficients();
		InitBcoefficients();
	}

	IIRLowPassFilter::~IIRLowPassFilter()
	{
	}

	void IIRLowPassFilter::ChangeCutoffFrequency(double newFpass)
	{
		assert(newFpass > 0);
		cutoff_frequency_ = newFpass;
		SetAngularCutoffFrequency();
		InitMemory();
		InitAcoefficients();
		InitBcoefficients();

	}

	double IIRLowPassFilter::GetCoefficientDenominator() const
	{
		return 4 * absorbtion_factor_ * absorbtion_factor_ * sample_rate_ * sample_rate_ +
			2 * absorbtion_factor_ * angular_cutoff_frequency_ * sample_rate_ * sqrt(double(2)) +
			angular_cutoff_frequency_ * angular_cutoff_frequency_;
	}

	double IIRLowPassFilter::GetA1CoefficientValue() const
	{
		auto denominator = GetCoefficientDenominator();
		return (2 * angular_cutoff_frequency_ * angular_cutoff_frequency_ -
			8 * absorbtion_factor_ * absorbtion_factor_ * sample_rate_ * sample_rate_) / denominator;
	}

	int IIRLowPassFilter::GetA2Coefficient() const
	{
		auto denominator = GetCoefficientDenominator();
		return (angular_cutoff_frequency_ * angular_cutoff_frequency_ +
			4 * absorbtion_factor_ * absorbtion_factor_ * sample_rate_ * sample_rate_ -
			2 * absorbtion_factor_ * angular_cutoff_frequency_ * sample_rate_ * sqrt(double(2))) / denominator;
	}

	float IIRLowPassFilter::FilterOutputLeft(float sample)
	{
		if (bypass_state_ == off)
			return FilterOutput(memory_left_, sample);
		else
			return sample;
	}

	float IIRLowPassFilter::FilterOutputRight(float sample)
	{
		if (bypass_state_ == off)
			return FilterOutput(memory_right_, sample);
		else
			return sample;
	}

	double IIRLowPassFilter::Spectrum(double frequency)
	{
		return sqrt((pow((b_coefficients_[0] + b_coefficients_[1] * cos(2 * PI*frequency) + b_coefficients_[2] * cos(4 * PI*frequency)), 2)
			+ pow(b_coefficients_[1] * sin(2 * PI*frequency) + b_coefficients_[2] * sin(4 * PI*frequency), 2)) /
			(pow((1 + a_coefficients_[1] * cos(2 * PI*frequency) + a_coefficients_[2] * cos(4 * PI*frequency)), 2)
				+ pow(a_coefficients_[1] * sin(2 * PI*frequency) + a_coefficients_[2] * sin(4 * PI*frequency), 2)));
	}

	void IIRLowPassFilter::InitAcoefficients()
	{
		if (a_coefficients_ == nullptr)
			a_coefficients_ = new double[COEFFICIENTS_NUMBER];



		a_coefficients_[0] = 1;
		a_coefficients_[1] = GetA1CoefficientValue();
		a_coefficients_[2] = GetA2Coefficient();
	}

	double IIRLowPassFilter::GetB0Coefficient() const
	{
		return (angular_cutoff_frequency_ * angular_cutoff_frequency_) / GetCoefficientDenominator();
	}

	void IIRLowPassFilter::SetAbsorbtionFactor(double absorbiton_db)
	{
		auto expression = pow(10, absorbiton_db / 10) - 1;
		absorbtion_factor_ = pow(expression, 1.0 / 4.0);
	}

	float IIRLowPassFilter::FilterOutput(float* memory, float sample) const
	{

		auto tmp = (sample - a_coefficients_[1] * memory[0] - a_coefficients_[2] * memory[1]);
		auto ans = tmp*b_coefficients_[0] + b_coefficients_[1] * memory[0] + b_coefficients_[2] * memory[1];
		memory[1] = memory[0];
		memory[0] = float(tmp);

		return float(ans);
	}

	void IIRLowPassFilter::InitBcoefficients()
	{
		if (b_coefficients_ == nullptr)
			b_coefficients_ = new double[MEMORY_SIZE];

		b_coefficients_[0] = GetB0Coefficient();
		b_coefficients_[1] = 2 * GetB0Coefficient();
		b_coefficients_[2] = GetB0Coefficient();
	}

	void IIRLowPassFilter::InitMemory()
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

	FirFilter::FirFilter(double sample_rate, double pass_band_frequency)
		: DigitalFilter(pass_band_frequency, sample_rate), filter_size_(0)

	{


	}

	FirFilter::FirFilter(double sample_rate, double pass_band_frequency, double* b_coefficients, int filter_size)
		: DigitalFilter(pass_band_frequency, sample_rate), filter_size_(filter_size), tmp_b_coeffcients_(b_coefficients)
	{
		InitMemory();
		InitBcoefficients();
	}

	FirFilter::~FirFilter()
	{
	}

	float FirFilter::FilterOutputLeft(float sample)
	{
		return FilterOutput(sample, memory_left_);
	}

	float FirFilter::FilterOutputRight(float sample)
	{
		return FilterOutput(sample, memory_right_);
	}

	double FirFilter::Spectrum(double frequency)
	{
		assert(b_coefficients_ != nullptr);
		auto real = 0.0;
		auto imag = 0.0;

		for (auto i = 0; i<filter_size_; ++i)
		{
			real += cos(2 * PI*frequency / sample_rate_)*b_coefficients_[i];
			imag += sin(2 * PI*frequency / sample_rate_)*b_coefficients_[i];
		}

		return sqrt(real*real + imag*imag);
	}

	void FirFilter::ChangeCutoffFrequency(double newFpass)
	{
	}

	void FirFilter::InitBcoefficients()
	{
		assert(filter_size_ != 0);
		if (b_coefficients_ != nullptr)
			delete b_coefficients_;
		b_coefficients_ = new double[filter_size_];

		for (int i = 0; i < filter_size_; ++i)
			b_coefficients_[i] = tmp_b_coeffcients_[i];

	}

	void FirFilter::InitAcoefficients()
	{
		a_coefficients_ = nullptr;
	}




	void FirFilter::ShiftMemory(float* memory) const
	{
		auto tmp = memory[0];
		for (auto i = 1; i<filter_size_; ++i)
		{
			memory[i] = memory[i - 1];
		}
	}

	void FirFilter::ClearMemory() const
	{
		for (auto i = 0; i<filter_size_; ++i)
		{
			memory_left_[i] = 0;
			memory_right_[i] = 0;
		}
	}

	void FirFilter::InitMemory()
	{
		assert(filter_size_ != 0);
		if (memory_left_ == nullptr &&
			memory_right_ == nullptr)
		{
			memory_left_ = new float[filter_size_];
			memory_right_ = new float[filter_size_];
		}
		ClearMemory();
	}

	float FirFilter::FilterOutput(float sample, float* memory) const
	{
		assert(memory != nullptr);
		assert(b_coefficients_ != nullptr);

		memory[0] = sample;
		double output = 0;

		for (auto i = 0; i < filter_size_; ++i)
			output += memory[i] * b_coefficients_[i];

		ShiftMemory(memory);

		return output;
	}

	void FirLowPassFilter::InitFilter()
	{
		InitDFactor();
		CalculateFilterSize();
		InitMemory();
		InitBetaFactor();
		InitBcoefficients();
	}

	FirLowPassFilter::FirLowPassFilter(double sample_rate, double pass_band_frequency, double stop_band_frequency, double absorbtion_in_stop_band)
		: FirFilter(sample_rate, pass_band_frequency), stop_band_frequency_(stop_band_frequency), absorbtion_in_stop_band_(absorbtion_in_stop_band)
	{
		InitFilter();

	}

	FirLowPassFilter::~FirLowPassFilter()
	{
	}

	void FirLowPassFilter::InitBcoefficients()
	{
		double h, w;

		if (b_coefficients_ != nullptr)
			delete b_coefficients_;
		b_coefficients_ = new double[filter_size_];
		auto fc = (cutoff_frequency_ / 2 + stop_band_frequency_ / 2) / sample_rate_;
		double half_of_filter_size = filter_size_ / 2;
		auto besseli_value = BesselZeroKindFunction(beta_factor_);
		for (int i = 0; i < filter_size_; ++i)
		{
			if (i - half_of_filter_size == 0)
				h = 2 * fc;
			else
				h = 2 * fc*sin(2 * 3.14*fc*(i - half_of_filter_size)) / (2 * 3.14*fc*(i - half_of_filter_size));
			w = BesselZeroKindFunction(beta_factor_*sqrt(1 - pow((i - half_of_filter_size) / half_of_filter_size, 2))) / besseli_value;
			b_coefficients_[i] = (w*h);
		}
	}

	void FirLowPassFilter::ChangeCutoffFrequency(double passband_requency, double stopband_frequency)
	{
		cutoff_frequency_ = passband_requency;
		stop_band_frequency_ = stopband_frequency;
		InitFilter();
	}

	double FirLowPassFilter::factorial(unsigned arg)
	{
		if (arg == 0 || arg == 1) return 1;
		int ans = 2;
		for (unsigned int i = 3; i <= arg; ++i)
			ans *= i;
		return ans;
	}

	void FirLowPassFilter::ChangeCutoffFrequency(double newFpass)
	{
		cutoff_frequency_ = newFpass;
		stop_band_frequency_ = 2 * cutoff_frequency_;
		InitFilter();
	}

	void FirLowPassFilter::CalculateFilterSize()
	{
		assert(d_factor_ != 0);
		assert(stop_band_frequency_ > cutoff_frequency_);
		auto df = stop_band_frequency_ - cutoff_frequency_;
		filter_size_ = static_cast<int>(ceil((d_factor_ * sample_rate_) / (df)));
		filter_size_ += (filter_size_ % 2 == 0) ? 1 : 0;
	}

	void FirLowPassFilter::InitDFactor()
	{
		assert(absorbtion_in_stop_band_ != 0);
		d_factor_ = (absorbtion_in_stop_band_ > 21) ? (absorbtion_in_stop_band_ - 7.95) / 14.36 : 0.922;
	}

	void FirLowPassFilter::InitBetaFactor()
	{
		assert(absorbtion_in_stop_band_ != 0);
		if (absorbtion_in_stop_band_ < 21)
		{
			beta_factor_ = 0.0;
			return;
		}
		if (absorbtion_in_stop_band_ >= 21 && absorbtion_in_stop_band_ <= 51)
		{
			beta_factor_ = 0.5842*pow(absorbtion_in_stop_band_ - 21, 0.4) + 0.07886*(absorbtion_in_stop_band_ - 21);
			return;
		}
		beta_factor_ = 0.1102*(absorbtion_in_stop_band_ - 8.7);
	}

	double FirLowPassFilter::BesselZeroKindFunction(double beta) const
	{
		double ans = 1;
		for (int k = 1; k < 7; ++k)
		{
			ans += pow(pow(beta / 2, k) / factorial(k), 2);
		}
		return ans;
	}



  void FirBandPassFilter::InitFilter()
  {
    InitDFactor();
    CalculateFilterSize();
    InitMemory();
    InitBetaFactor();
    InitBcoefficients();
  }

  FirBandPassFilter::FirBandPassFilter(double sample_rate, double pass_band_frequency, double stop_band_frequency, double absorbtion_in_stop_band, double middle_frequency)
    : FirFilter(sample_rate, pass_band_frequency), stop_band_frequency_(stop_band_frequency), absorbtion_in_stop_band_(absorbtion_in_stop_band), middle_frequency_(middle_frequency)
  {
    InitFilter();

  }

  FirBandPassFilter::~FirBandPassFilter()
  {

  }


  void FirBandPassFilter::InitBcoefficients()
  {
    double h, w, f0;

    if (b_coefficients_ != nullptr)
      delete b_coefficients_;
    b_coefficients_ = new double[filter_size_];
    auto fc = (cutoff_frequency_ / 2 + stop_band_frequency_ / 2) / sample_rate_;
    double half_of_filter_size = filter_size_ / 2;
    auto besseli_value = BesselZeroKindFunction(beta_factor_);
    for (int i = 0; i < filter_size_; ++i)
    {
      if (i - half_of_filter_size == 0)
        h = 2 * fc;
      else
        h = 2 * fc*sin(2 * 3.14*fc*(i - half_of_filter_size)) / (2 * 3.14*fc*(i - half_of_filter_size));
      w = BesselZeroKindFunction(beta_factor_*sqrt(1 - pow((i - half_of_filter_size) / half_of_filter_size, 2))) / besseli_value;
      f0 = std::cos(2 * PI*middle_frequency_*i / sample_rate_);
      b_coefficients_[i] = (w*h*f0);
    }
  }

  void FirBandPassFilter::ChangeCutoffFrequency(double passband_requency, double stopband_frequency, double middle_frequency)
  {
    cutoff_frequency_ = passband_requency;
    stop_band_frequency_ = stopband_frequency;
    middle_frequency_ = middle_frequency;
    InitFilter();
  }

  double FirBandPassFilter::factorial(unsigned arg)
  {
    if (arg == 0 || arg == 1) return 1;
    int ans = 2;
    for (unsigned int i = 3; i <= arg; ++i)
      ans *= i;
    return ans;
  }

  void FirBandPassFilter::ChangeCutoffFrequency(double newFpass)
  {
    cutoff_frequency_ = newFpass;
    stop_band_frequency_ = 2 * cutoff_frequency_;
    InitFilter();
  }

  void FirBandPassFilter::CalculateFilterSize()
  {
    assert(d_factor_ != 0);
    assert(stop_band_frequency_ > cutoff_frequency_);
    auto df = stop_band_frequency_ - cutoff_frequency_;
    filter_size_ = static_cast<int>(ceil((d_factor_ * sample_rate_) / (df)));
    filter_size_ += (filter_size_ % 2 == 0) ? 1 : 0;
  }

  void FirBandPassFilter::InitDFactor()
  {
    assert(absorbtion_in_stop_band_ != 0);
    d_factor_ = (absorbtion_in_stop_band_ > 21) ? (absorbtion_in_stop_band_ - 7.95) / 14.36 : 0.922;
  }

  void FirBandPassFilter::InitBetaFactor()
  {
    assert(absorbtion_in_stop_band_ != 0);
    if (absorbtion_in_stop_band_ < 21)
    {
      beta_factor_ = 0.0;
      return;
    }
    if (absorbtion_in_stop_band_ >= 21 && absorbtion_in_stop_band_ <= 51)
    {
      beta_factor_ = 0.5842*pow(absorbtion_in_stop_band_ - 21, 0.4) + 0.07886*(absorbtion_in_stop_band_ - 21);
      return;
    }
    beta_factor_ = 0.1102*(absorbtion_in_stop_band_ - 8.7);
  }

  double FirBandPassFilter::BesselZeroKindFunction(double beta) const
  {
    double ans = 1;
    for (int k = 1; k < 7; ++k)
    {
      ans += pow(pow(beta / 2, k) / factorial(k), 2);
    }
    return ans;
  }



	FirPolyphaseDecimatorFilter::FirPolyphaseDecimatorFilter(double sample_rate, int decimation_factor)
		:FirLowPassFilter(sample_rate, sample_rate / (3 * decimation_factor), sample_rate / decimation_factor, 80), decimation_factor_(decimation_factor), number_of_polyphase_filters(decimation_factor),
		sample_counter_left_(-1), sample_counter_right_(-1)
	{
		InitDecimatorFilter();
	}

	FirPolyphaseDecimatorFilter::~FirPolyphaseDecimatorFilter()
	{
	}

	void FirPolyphaseDecimatorFilter::CreatePolyphaseFilter(double* tmp_b_coefficients, int i)
	{
		for (int j = 0; j < polyphase_filter_size_; ++j)
		{
			if (i + j*decimation_factor_ < filter_size_)
				tmp_b_coefficients[j] = b_coefficients_[i + j*decimation_factor_];
			else
				tmp_b_coefficients[j] = 0;
		}
		polyphase_filters_.push_back(new FirFilter(sample_rate_ / decimation_factor_, 1, tmp_b_coefficients, polyphase_filter_size_));
	}

	void FirPolyphaseDecimatorFilter::InitPolyphaseFilters()
	{
		double* tmp_b_coefficients = new double[polyphase_filter_size_];
		for (int i = 0; i < number_of_polyphase_filters; ++i)
		{

			CreatePolyphaseFilter(tmp_b_coefficients, i);
		}
		delete[] tmp_b_coefficients;
	}

	void FirPolyphaseDecimatorFilter::InitFilterOutpuTables()
	{
		output_left_ = new float[number_of_polyphase_filters];
		output_right_ = new float[number_of_polyphase_filters];

		for (int i = 0; i<number_of_polyphase_filters; ++i)
		{
			output_left_[i] = 0;
			output_right_[i] = 0;
		}
	}

	void FirPolyphaseDecimatorFilter::InitDecimatorFilter()
	{

		assert(filter_size_ != 0);


		polyphase_filter_size_ = int(ceil(double(filter_size_) / decimation_factor_));

		InitFilterOutpuTables();


		InitPolyphaseFilters();
	}

	float FirPolyphaseDecimatorFilter::FilterOutputLeft(float sample)
	{
		++sample_counter_left_;
		sample_counter_left_ %= decimation_factor_;
		if (sample_counter_left_ == 0)
			output_right_[sample_counter_left_] = polyphase_filters_[sample_counter_left_]->FilterOutputLeft(sample);
		else
		{
			int index = number_of_polyphase_filters - sample_counter_left_;
			output_left_[index] = polyphase_filters_[index]->FilterOutputLeft(sample);
		}

		float ans = 0;
		for (int i = 0; i < number_of_polyphase_filters; ++i)
			ans += output_left_[i];

		return ans;

	}

	float FirPolyphaseDecimatorFilter::FilterOutputRight(float sample)
	{
		++sample_counter_right_;
		sample_counter_right_ %= decimation_factor_;
		if (sample_counter_right_ == 0)
			output_right_[sample_counter_right_] = polyphase_filters_[sample_counter_right_]->FilterOutputRight(sample);
		else
		{
			int index = number_of_polyphase_filters - sample_counter_right_ - 1;
			output_right_[index] = polyphase_filters_[index]->FilterOutputRight(sample);
		}

		float ans = 0;
		for (int i = 0; i < number_of_polyphase_filters; ++i)
			ans += output_right_[i];

		return ans;
	}

	float FirPolyphaseDecimatorFilter::FilterOutput(float sample, int& sample_counter)
	{
		if (sample_counter == 0)
			return polyphase_filters_[sample_counter]->FilterOutputLeft(sample);
		else
			return polyphase_filters_[polyphase_filter_size_ - sample_counter - 1]->FilterOutputLeft(sample);
	}


	IIRBandPassFilter::IIRBandPassFilter(double cutoff_frequency, double sample_rate, double gain, double q_factor, bypassState bypass_state)
		: IIRParametricBandPassFilter(cutoff_frequency, sample_rate, gain, q_factor, bypass_state)
	{
    ChangeQFactor(q_factor);
	}

	float IIRBandPassFilter::FilterOutputLeft(float sample)
	{
		if (bypass_state_ != off)
			return sample;
		double ans;
		auto tmp = (sample - a_coefficients_[1] * memory_left_[0] - a_coefficients_[2] * memory_left_[1]);
		ans = tmp*b_coefficients_[0] + b_coefficients_[1] * memory_left_[0] + b_coefficients_[2] * memory_left_[1];
		memory_left_[1] = memory_left_[0];
		memory_left_[0] = float(tmp);


		return float(ans);
	}

	float IIRBandPassFilter::FilterOutputRight(float sample)
	{
		if (bypass_state_ != off)
			return sample;

		double ans;
		auto tmp = (sample - a_coefficients_[1] * memory_right_[0] - a_coefficients_[2] * memory_right_[1]);
		ans = tmp*b_coefficients_[0] + b_coefficients_[1] * memory_right_[0] + b_coefficients_[2] * memory_right_[1];
		memory_right_[1] = memory_right_[0];
		memory_right_[0] = float(tmp);


		return float(ans);
	}

	double IIRBandPassFilter::Spectrum(double frequency)
	{
		if (bypass_state_ == off)
		{
			double b0 = b_coefficients_[0], b1 = b_coefficients_[1], b2 = b_coefficients_[2];

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

	void IIRBandPassFilter::InitBcoefficients()
	{
		if (b_coefficients_ == nullptr)
			b_coefficients_ = new double[COEFFICIENTS_NUMBER];

		auto denominator = 4 * sample_rate_*sample_rate_ + angular_cutoff_frequency_*angular_cutoff_frequency_ + 2 * bandwidth_*sample_rate_;
		b_coefficients_[0] = 2 * bandwidth_*(sample_rate_ / denominator);
		b_coefficients_[1] = 0;
		b_coefficients_[2] = -2 * bandwidth_*(sample_rate_ / denominator);
	}
}