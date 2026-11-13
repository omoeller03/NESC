#pragma once

#include "nesc/Nucleus.hpp"

namespace nesc
{

	enum class ChargeDistribution
	{
		UniformSphere,
		Fermi
	};

	class FormFactor
	{
		public:
			explicit FormFactor(
					const Nucleus& nucleus,
					ChargeDistribution model = ChargeDistribution::UniformSphere,
					double diffuseness_fm = 0.5
					);

			// Calculate the form factor for a given momentum transfer q [fm^-1].
			double evaluate(double q_fm_inv) const;

			// Analytic form factor of a uniformly charged sphere.
			double uniform_sphere(double q_fm_inv) const;

			// Numerically integrated Fermi charge distribution.
			double fermi(double q_fm_inv) const;

			// Parameters of the selected charge distribution.
			double radius_parameter_fm() const;
			double diffuseness_fm() const;

		private:
			double rms_charge_radius_fm_;
			double radius_parameter_fm_;
			double diffuseness_fm_;

			ChargeDistribution model_;

			/*
			double fermi_integral(
					double q_fm_inv,
					double power,
					double upper_limit_fm
					) const;
			*/

			double fermi_integral(
					double q_fm_inv,
					double power,
					double c_fm,
					double upper_limit_fm
					) const;

			double fermi_rms_radius_squared(double c_fm) const;

			double determine_fermi_radius_parameter() const;
	};

} // namespace nesc



/*
#pragma once

#include "nesc/Nucleus.hpp"

namespace nesc
{

	class FormFactor {
		public:
			explicit FormFactor(const Nucleus& nucleus);

			double uniform_sphere(double q_fm_inv) const;

		private:
			double radius_fm_;
	};

} // namespace nesc
*/
