#include "nesc/Constants.hpp"
#include "nesc/CrossSection.hpp"
#include "nesc/Kinematics.hpp"
#include "nesc/Nucleus.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{

	struct Options
	{
		int Z = 0;
		int A = 0;
		double energy_MeV = 0.0;

		bool radius_provided = false;
		double radius_fm = 0.0;

		double theta_start_deg = 10.0;
		double theta_end_deg = 150.0;
       		double theta_step_deg = 1.0;

		bool plot = false;
	};

	void print_help(const char* program)
	{
		std::cout
			<< "NESC - Nuclear Electron Scattering Calculator\n\n"
			<< "Usage:\n"
			<< "  " << program
			<< " --Z Z --A A --energy MeV --radius fm"
			<< " [--theta start:end:step] [--plot]\n\n"
			<< "Options:\n"
			<< "  --Z       Nuclear charge number\n"
			<< "  --A       Nuclear mass number\n"
			<< "  --energy  Incident electron kinetic energy in MeV\n"
			//<< "  --radius  RMS nuclear charge radius in fm\n"
			<< "  --radius  RMS nuclear charge radius in fm (optional; estimated from A if omitted)\n"
			<< "  --theta   Angle range in degrees (default: 10:150:1)\n"
			<< "  --plot    Open the result with Gnuplot\n"
			<< "  --help    Show this message\n\n"
			<< "Model: elastic PWBA, uniform spherical charge distribution.\n";
	}

	double parse_double(const std::string& value, const std::string& name)
	{
		std::size_t used = 0;
		const double result = std::stod(value, &used);

		if(used != value.size() || !std::isfinite(result))
		{
			throw std::invalid_argument("Invalid value for " + name);
		}

		return result;
	}

	int parse_int(const std::string& value, const std::string& name)
	{
		std::size_t used = 0;
		const int result = std::stoi(value, &used);

		if(used != value.size())
		{
			throw std::invalid_argument("Invalid value for " + name);
		}

		return result;
	}

	void parse_theta_range(const std::string& value, Options& options)
	{
		std::istringstream input(value);
		std::string start, end, step;

		if(!std::getline(input, start, ':') ||
				!std::getline(input, end, ':') ||
				!std::getline(input, step, ':')) {
			throw std::invalid_argument("Theta range must have the form start:end:step.");
		}

		std::string extra;
		if(std::getline(input, extra, ':'))
		{
			throw std::invalid_argument("Theta range must have exactly three values.");
		}

		options.theta_start_deg = parse_double(start, "--theta");
		options.theta_end_deg   = parse_double(end, "--theta");
		options.theta_step_deg  = parse_double(step, "--theta");

		if(options.theta_start_deg <= 0.0 ||
				options.theta_end_deg > 180.0 ||
				options.theta_end_deg < options.theta_start_deg ||
				options.theta_step_deg <= 0.0) {
			throw std::invalid_argument("Require 0 < start <= end <= 180 and step > 0.");
		}
	}

	Options parse_arguments(int argc, char* argv[])
	{
		Options options;

		for(int i = 1; i < argc; ++i)
		{
			const std::string arg = argv[i];

			if(arg == "--help" || arg == "-h") {
				print_help(argv[0]);
				std::exit(0);
			} else if (arg == "--plot") {
				options.plot = true;
			} else if (arg == "--Z" || arg == "--A" || arg == "--energy" || arg == "--radius" || arg == "--theta") {
				if(i + 1 >= argc)
				{
					throw std::invalid_argument("Missing value after " + arg);
				}

				const std::string value = argv[++i];

				if(arg == "--Z") {
					options.Z = parse_int(value, arg);
				} else if (arg == "--A") {
					options.A = parse_int(value, arg);
				} else if (arg == "--energy") {
					options.energy_MeV = parse_double(value, arg);
				} else if (arg == "--radius") {
					options.radius_fm = parse_double(value, arg);
					options.radius_provided = true;
				} else if (arg == "--theta") {
					parse_theta_range(value, options);
				}
			} else {
				throw std::invalid_argument("Unknown option: " + arg);
			}
		}

		if(options.Z == 0 || options.A == 0 || options.energy_MeV <= 0.0)
		{
			throw std::invalid_argument("Please specify --Z, --A, --energy and optionally --radius.");
		}

		if(options.radius_provided && options.radius_fm <= 0.0)
		{
			throw std::invalid_argument("--radius must be positive if specified.");
		}

		return options;

	}

} // namespace

int main(int argc, char* argv[])
{

	try
	{
		const Options options = parse_arguments(argc, argv);



		double radius_fm = options.radius_fm;
		if(!options.radius_provided)
		{
			constexpr double r0_fm = 1.2;

			radius_fm = std::sqrt(3.0 / 5.0)
				* r0_fm
				* std::cbrt(static_cast<double>(options.A));

			std::cout
				<< "No charge radius specified.\n"
				<< "Using estimated RMS charge radius: "
				<< radius_fm << " fm\n";
		} else {
			std::cout
				<< "Using user-specified RMS charge radius: "
				<< radius_fm << " fm\n";
		}


		nesc::Nucleus nucleus
		{
			options.Z,
			options.A,
			radius_fm
		};

		nesc::Kinematics kinematics
		{
			options.energy_MeV
		};

		nucleus.validate();
		kinematics.validate();

		const nesc::CrossSection cross_section(nucleus, kinematics);

		std::ofstream output("cross_section.dat");
		if(!output)
		{
			throw std::runtime_error("Could not open cross_section.dat for writing.");
		}

		output << "# NESC: elastic PWBA, uniform spherical charge distribution\n";
		output << "# Z = " << nucleus.Z
			<< ", A = " << nucleus.A
			<< ", T_e = " << kinematics.kinetic_energy_MeV << " MeV"
			<< ", rms charge radius = "
			<< nucleus.rms_charge_radius_fm << " fm\n";
		output << "# theta_deg  q_fm^-1  F(q)  dsigma_dOmega_mb_sr\n";
		output << std::scientific << std::setprecision(10);

		const double span = options.theta_end_deg - options.theta_start_deg;
		const int number_of_steps = static_cast<int>(std::floor(span / options.theta_step_deg + 1.0e-10));

		for(int i = 0; i <= number_of_steps; ++i)
		{
			const double theta_deg = options.theta_start_deg + i * options.theta_step_deg;

			const double theta_rad = theta_deg * nesc::constants::pi / 180.0;

			output << theta_deg << " "
				<< cross_section.momentum_transfer(theta_rad) << " "
				<< cross_section.form_factor(theta_rad) << " "
				<< cross_section.pwba_mb_sr(theta_rad) << "\n";
		}

		output.close();

		std::cout << "Calculation complete.\n"
			<< "Output: cross_section.dat\n"
			<< "Model: elastic PWBA, uniform spherical charge distribution\n";

		if(options.plot)
		{
			const int status = std::system("gnuplot -persist scripts/plot_elastic_cross_section.gp");

			if(status != 0)
			{
				std::cerr
					<< "Warning: Gnuplot could not be started.\n"
					<< "Check that Gnuplot is installed and available in PATH.\n";
				return 1;
			}
		}

		return 0;
	} 
	catch(const std::exception& error)
	{
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
