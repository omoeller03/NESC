#pragma once

namespace nesc
{

	struct Nucleus {
		int Z = 0;                    // Proton number
		int A = 0;                    // Mass number
		double rms_charge_radius_fm = 0.0;

		void validate() const;
	};

} // namespace nesc
