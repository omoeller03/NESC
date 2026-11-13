#include "nesc/FormFactor.hpp"
#include "nesc/Kinematics.hpp"
#include "nesc/Nucleus.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

} // namespace

int main()
{
    try {
        const nesc::Nucleus cerium{58, 140, 4.9};
        const nesc::FormFactor form_factor(cerium);

        require(
            std::abs(form_factor.uniform_sphere(0.0) - 1.0) < 1.0e-12,
            "F(0) must equal 1."
        );

        const double small_q = form_factor.uniform_sphere(1.0e-5);
        require(
            std::abs(small_q - 1.0) < 1.0e-7,
            "F(q) should approach 1 for small q."
        );

        const nesc::Kinematics kinematics{65.0};

        require(
            kinematics.momentum_MeV_c() > 65.0,
            "Electron momentum should be positive and relativistic."
        );

        require(
            kinematics.beta() > 0.0 && kinematics.beta() < 1.0,
            "Electron beta must lie between 0 and 1."
        );

        require(
            std::abs(kinematics.momentum_transfer_fm_inv(0.0)) < 1.0e-14,
            "Momentum transfer must vanish at zero angle."
        );

        std::cout << "All physics tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Test failed: " << error.what() << "\n";
        return 1;
    }
}
