#pragma once
#ifndef DIGITAL_FILTERS_H_
#define DIGITALFILTERS_H_

#define PI 3.141592653589793
#include <cassert>
#include <vector>
#include <assert.h>
#include <cmath>
namespace HoffFilters
{
	enum bypassState { on, off };



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


	

	class FIRInterpolatorFilter :
		public DigitalFilter
	{
	public:
		FIRInterpolatorFilter(int filter_size, int interpolation_factor);
		~FIRInterpolatorFilter();

		float FilterOutput(float * samples, int start, int samples_vector_size) const;

		void InitBcoefficients() override;
		void InitMemory() override;
		void ShiftMemory(float sample) const;

		void ChangeCutoffFrequency(double newFpass) override {};
		float FilterOutputLeft(float sample) override { return 0; };
		float FilterOutputRight(float sample) override { return 0; };
		double Spectrum(double frequency) override { return 0; };
		void InitAcoefficients() override {};

	protected:

		int filter_size_;
		int interpolation_factor_;
	};

	


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
		float FilterOutput(float* memory, float sample) const;

		static const int MEMORY_SIZE = 2, COEFFICIENTS_NUMBER = 3;
		double absorbtion_factor_;
	};




	class FirFilter :
		public DigitalFilter
	{
	public:
		FirFilter(double sample_rate, double pass_band_frequency);
		FirFilter(double sample_rate, double pass_band_frequency, double* b_coefficients, int filter_size);
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


	



	class FirLowPassFilter :
		public FirFilter
	{
	public:
		void InitFilter();
		FirLowPassFilter(double sample_rate, double pass_band_frequency, double stop_band_frequency, double absorbtion_in_stop_band);
		~FirLowPassFilter();

		virtual void InitBcoefficients() override;
		void ChangeCutoffFrequency(double passband_requency, double stopband_frequency);

		void CalculateFilterSize();
		void InitDFactor();
		void InitBetaFactor();
    virtual float FilterOutputLeft(float sample) override;
    virtual float FilterOutputRight(float sample) override;
	protected:
    float FilterOutput(float sample, float* memory) const;
		double stop_band_frequency_;
		double absorbtion_in_stop_band_;
		double beta_factor_;
		double d_factor_;
	private:
		static double factorial(unsigned int arg);
		void ChangeCutoffFrequency(double newFpass) override;
		double BesselZeroKindFunction(double beta) const;


	};

  class FirBandPassFilter :
    public FirFilter
  {
  public:
    void InitFilter();
    FirBandPassFilter(double sample_rate, double pass_band_down, double stop_band_down,double pass_band_up, double stop_band_up,  double absorbtion_in_stop_band);
    ~FirBandPassFilter();

    virtual void InitBcoefficients() override;
    void ChangeCutoffFrequency(double pass_band_down, double stop_band_down, double pass_band_up, double stop_band_up);

    void CalculateFilterSize();
    void InitDFactor();
    void InitBetaFactor();
  protected:
    double stop_band_down_;
    double absorbtion_in_stop_band_;
    double pass_band_up_;
    double stop_band_up_;
    double delta_frequency_;
    double beta_factor_;
    double d_factor_;
  private:
    static double factorial(unsigned int arg);
    void SetDeltaFrequency();
    void ChangeCutoffFrequency(double newFpass) override;
    double BesselZeroKindFunction(double beta) const;
  };


	class FirPolyphaseDecimatorFilter : FirLowPassFilter
	{
	public:
		FirPolyphaseDecimatorFilter(double sample_rate, int decimation_factor);
		~FirPolyphaseDecimatorFilter();
		void CreatePolyphaseFilter(double* tmp_b_coefficients, int i);
		void InitPolyphaseFilters();
		void InitFilterOutpuTables();
		void InitDecimatorFilter();

		virtual float FilterOutputLeft(float sample) override;
		virtual float FilterOutputRight(float sample) override;

	private:
		float FilterOutput(float sample, int& sample_counter);
		float* output_left_;
		float* output_right_;
		const int decimation_factor_;
		std::vector<FirFilter*> polyphase_filters_;
		int number_of_polyphase_filters;
		int polyphase_filter_size_;
		int sample_counter_left_;
		int sample_counter_right_;

	};

	
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

}

#endif
