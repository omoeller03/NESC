#include "nesc/CrossSection.hpp"
#include "nesc/Constants.hpp"

#include <cmath>
#include <stdexcept>

namespace nesc {

CrossSection::CrossSection(
    const Nucleus& nucleus,
    const Kinematics& kinematics)
    : nucleus_(nucleus),
      kinematics_(kinematics),
      form_factor_(nucleus)
{
    nucleus_.validate();
    kinematics_.validate();
}

double CrossSection::mott_fm2_sr(double theta_rad) const
{
    if (!std::isfinite(theta_rad) ||
        theta_rad <= 0.0 ||
        theta_rad > constants::pi) {
        throw std::invalid_argument(
            "Mott cross section requires 0 < theta <= pi."
        );
    }

    const double half = theta_rad / 2.0;
    const double sin_half = std::sin(half);
    const double beta = kinematics_.beta();
    const double energy = kinematics_.total_energy_MeV();

    const double numerator =
        nucleus_.Z * constants::alpha
        * constants::hbar_c_MeV_fm;

    const double denominator =
        2.0 * energy * beta * beta
        * sin_half * sin_half;

    const double spin_factor =
        1.0 - beta * beta * sin_half * sin_half;

    return (numerator / denominator)
         * (numerator / denominator)
         * spin_factor;
}

double CrossSection::form_factor(double theta_rad) const
{
    return form_factor_.uniform_sphere(
        momentum_transfer(theta_rad)
    );
}

double CrossSection::momentum_transfer(double theta_rad) const
{
    return kinematics_.momentum_transfer_fm_inv(theta_rad);
}

double CrossSection::pwba_fm2_sr(double theta_rad) const
{
    const double F = form_factor(theta_rad);
    return mott_fm2_sr(theta_rad) * F * F;
}

double CrossSection::pwba_mb_sr(double theta_rad) const
{
    return pwba_fm2_sr(theta_rad) * constants::fm2_to_mb;
}

} // namespace nesc
