#pragma once

#include "nesc/Kinematics.hpp"
#include "nesc/Nucleus.hpp"
#include "nesc/FormFactor.hpp"

namespace nesc
{

	class CrossSection
	{
		public:
			CrossSection(const Nucleus& nucleus, const Kinematics& kinematics);

			double mott_fm2_sr(double theta_rad) const;

			double pwba_fm2_sr(double theta_rad) const;

			double pwba_mb_sr(double theta_rad) const;

			double form_factor(double theta_rad) const;

			double momentum_transfer(double theta_rad) const;

		private:
			Nucleus nucleus_;
			Kinematics kinematics_;
			FormFactor form_factor_;
};

} // namespace nesc
