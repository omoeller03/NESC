#pragma once

namespace nesc
{

	struct Kinematics
	{
		double kinetic_energy_MeV = 0.0;

		double total_energy_MeV() const;
		double momentum_MeV_c() const;
		double beta() const;

		double momentum_transfer_fm_inv(double theta_rad) const;
		void validate() const;
	};

} // namespace nesc
