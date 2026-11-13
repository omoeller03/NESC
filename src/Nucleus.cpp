#include "nesc/Nucleus.hpp"

#include <stdexcept>

namespace nesc
{

	void Nucleus::validate() const
	{
		if(Z <= 0)
		{
			throw std::invalid_argument("Z must be greater than zero.");
		}
	
		if(A < Z)
		{
			throw std::invalid_argument("A must be greater than or equal to Z.");
		}
	
		if(rms_charge_radius_fm <= 0.0)
		{
			throw std::invalid_argument("The RMS charge radius must be greater than zero.");
		}
	}

} // namespace nesc
