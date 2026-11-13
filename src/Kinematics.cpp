#include "nesc/Kinematics.hpp"
#include "nesc/Constants.hpp"

#include <cmath>
#include <stdexcept>

namespace nesc
{

	void Kinematics::validate() const
	{
		if (!std::isfinite(kinetic_energy_MeV) || kinetic_energy_MeV <= 0.0) {
			throw std::invalid_argument("Electron kinetic energy must be finite and positive.");
		}
	}

	double Kinematics::total_energy_MeV() const
	{
		validate();
		return kinetic_energy_MeV + constants::electron_mass_MeV;
	}

	double Kinematics::momentum_MeV_c() const
	{
		validate();

		const double m = constants::electron_mass_MeV;
		const double T = kinetic_energy_MeV;

		// (pc)^2 = T(T + 2 mc^2), with c = 1 in MeV units.
		return std::sqrt(T * (T + 2.0 * m));
	}

	double Kinematics::beta() const
	{
		return momentum_MeV_c() / total_energy_MeV();
	}

	double Kinematics::momentum_transfer_fm_inv(double theta_rad) const
	{
		if(!std::isfinite(theta_rad) || theta_rad < 0.0 || theta_rad > constants::pi)
		{
			throw std::invalid_argument("Scattering angle must be between 0 and pi radians.");
		}

		return 2.0 * momentum_MeV_c() * std::sin(theta_rad / 2.0) / constants::hbar_c_MeV_fm;
	}

} // namespace nesc
