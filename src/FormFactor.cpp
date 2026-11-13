#include "nesc/FormFactor.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nesc
{

	namespace
	{

		constexpr double integration_tail = 40.0;
		constexpr int integration_intervals = 4000;
		constexpr int bisection_iterations = 100;

		double fermi_density_shape(double r_fm,double c_fm,double diffuseness_fm)
		{
			const double x = (r_fm - c_fm) / diffuseness_fm;

			// Avoid numerical overflow in exp(x).
			if (x > 700.0)
			{
				return 0.0;
			}

			if (x < -700.0)
			{
				return 1.0;
			}

			return 1.0 / (1.0 + std::exp(x));
		}

		double sinc(double x)
		{
			if(std::abs(x) < 1.0e-4)
			{
				const double x2 = x * x;
				return 1.0 - x2 / 6.0 + x2 * x2 / 120.0;
			}

			return std::sin(x) / x;
		}

	} // namespace




	FormFactor::FormFactor(const Nucleus& nucleus,ChargeDistribution model,double diffuseness_fm)
		: rms_charge_radius_fm_(nucleus.rms_charge_radius_fm),
		radius_parameter_fm_(0.0),
		diffuseness_fm_(diffuseness_fm),
		model_(model)
		
	{
		nucleus.validate();

		if(!std::isfinite(diffuseness_fm_) || diffuseness_fm_ <= 0.0)
		{
			throw std::invalid_argument("Diffuseness must be finite and positive.");
		}

		if(model_ == ChargeDistribution::UniformSphere) {
			radius_parameter_fm_ = std::sqrt(5.0 / 3.0) * rms_charge_radius_fm_;
		} else {
			radius_parameter_fm_ = determine_fermi_radius_parameter();
		}
	}

	double FormFactor::evaluate(double q_fm_inv) const
	{
		if(model_ == ChargeDistribution::UniformSphere)
		{
			return uniform_sphere(q_fm_inv);
		}

		return fermi(q_fm_inv);
	}

	double FormFactor::uniform_sphere(double q_fm_inv) const
	{
		if(!std::isfinite(q_fm_inv) || q_fm_inv < 0.0)
		{
			throw std::invalid_argument("Momentum transfer must be finite and non-negative.");
		}

		const double x = q_fm_inv * radius_parameter_fm_;

		if(std::abs(x) < 1.0e-4)
		{
			const double x2 = x * x;
			return 1.0 - x2 / 10.0 + x2 * x2 / 280.0;
		}

		return 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
	}

	double FormFactor::fermi(double q_fm_inv) const
	{
		if(!std::isfinite(q_fm_inv) || q_fm_inv < 0.0) {
			throw std::invalid_argument("Momentum transfer must be finite and non-negative.");
		}

		const double upper_limit = radius_parameter_fm_ + integration_tail * diffuseness_fm_;

		const double numerator = fermi_integral(q_fm_inv, 2.0, radius_parameter_fm_, upper_limit); //fermi_integral(q_fm_inv, 2.0, upper_limit);

		const double denominator = fermi_integral(0.0, 2.0, radius_parameter_fm_, upper_limit); //fermi_integral(0.0, 2.0, upper_limit);

		if(denominator <= 0.0)
		{
			throw std::runtime_error("Invalid normalization integral for Fermi distribution.");
		}

		return numerator / denominator;
	}

	/*
	double FormFactor::fermi_integral(double q_fm_inv,double power,double upper_limit_fm) const
	{
		// Composite Simpson integration.
		// The number of intervals must be even.
		constexpr int n = integration_intervals;

		const double h = upper_limit_fm / n;
		double sum = 0.0;

		for(int i = 0; i <= n; ++i)
		{
			const double r = i * h;

			const double density = fermi_density_shape(r,radius_parameter_fm_,diffuseness_fm_);

			double value = std::pow(r, power) * density;

			if(q_fm_inv != 0.0)
			{
				value *= sinc(q_fm_inv * r);
			}

			int weight;

			if(i == 0 || i == n) {
				weight = 1;
			} else if (i % 2 == 0) {
				weight = 2;
			} else {
				weight = 4;
			}

			sum += weight * value;
		}

		return sum * h / 3.0;
	}
	*/


	double FormFactor::fermi_integral(
    double q_fm_inv,
    double power,
    double c_fm,
    double upper_limit_fm) const
{
    constexpr int n = integration_intervals;
    const double h = upper_limit_fm / n;
    double sum = 0.0;

    for (int i = 0; i <= n; ++i) {
        const double r = i * h;

        const double density = fermi_density_shape(
            r, c_fm, diffuseness_fm_
        );

        double value = std::pow(r, power) * density;

        if (q_fm_inv != 0.0) {
            value *= sinc(q_fm_inv * r);
        }

        const int weight =
            (i == 0 || i == n) ? 1 :
            (i % 2 == 0) ? 2 : 4;

        sum += weight * value;
    }

    return sum * h / 3.0;
}




	double FormFactor::fermi_rms_radius_squared(double c_fm) const
	{
		const double upper_limit = c_fm + integration_tail * diffuseness_fm_;

		// The integration helper uses radius_parameter_fm_.
		// Save and restore it so the trial value c_fm is used here.
		const double original_c = radius_parameter_fm_;
		const_cast<FormFactor*>(this)->radius_parameter_fm_ = c_fm;

		const double integral_r2 = fermi_integral(0.0, 2.0, upper_limit);

		const double integral_r4 = fermi_integral(0.0, 4.0, upper_limit);

		const_cast<FormFactor*>(this)->radius_parameter_fm_ = original_c;

		if(integral_r2 <= 0.0)
		{
			throw std::runtime_error("Invalid Fermi RMS-radius integral.");
		}

		return integral_r4 / integral_r2;
	}


	double FormFactor::determine_fermi_radius_parameter() const
	{
		const double target = rms_charge_radius_fm_ * rms_charge_radius_fm_;

		double low = 0.0;
		double high = std::max(
				rms_charge_radius_fm_,
				4.0 * diffuseness_fm_
				);

		// Find an upper bracket for the target RMS radius.
		for(int i = 0; i < 100; ++i)
		{
			if(fermi_rms_radius_squared(high) >= target) { break; }

			high *= 2.0;
		}

		if(fermi_rms_radius_squared(high) < target)
		{
			throw std::runtime_error("Could not determine Fermi radius parameter.");
		}


		for(int i = 0; i < bisection_iterations; ++i) {
			const double mid = 0.5 * (low + high);

			if(fermi_rms_radius_squared(mid) < target) {
				low = mid;
			} else {
				high = mid;
			}
		}

		return 0.5 * (low + high);
	}

	double FormFactor::radius_parameter_fm() const
	{
		return radius_parameter_fm_;
	}

	double FormFactor::diffuseness_fm() const
	{
		return diffuseness_fm_;
	}

} // namespace nesc




/*
#include "nesc/FormFactor.hpp"

#include <cmath>
#include <stdexcept>

namespace nesc
{

	FormFactor::FormFactor(const Nucleus& nucleus)
	{
		nucleus.validate();

		radius_fm_ = std::sqrt(5.0 / 3.0) * nucleus.rms_charge_radius_fm;
	}

	double FormFactor::uniform_sphere(double q_fm_inv) const
	{
		if(!std::isfinite(q_fm_inv) || q_fm_inv < 0.0)
		{
			throw std::invalid_argument("Momentum transfer must be finite and non-negative.");
		}

		const double x = q_fm_inv * radius_fm_;

		// Avoid numerical cancellation in the q -> 0 limit.
		if(std::abs(x) < 1.0e-4)
		{
			const double x2 = x * x;
			return 1.0 - x2 / 10.0 + x2 * x2 / 280.0;
		}

		return 3.0 * (std::sin(x) - x * std::cos(x)) / (x * x * x);
	}

} // namespace nesc
*/
