#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

enum integration_method {    //use this type to specify the integration method
    euler,                  //simple euler step
    leap_frog               //better use leap frog
};

inline std::ostream &operator << (std::ostream &stream, const integration_method method)
{
    if (method == euler)
        stream << "eulerstep" << std::endl;
    else if (method == leap_frog)
        stream << "leapfrog integration" << std::endl;
    else {
        std::cerr << "unkown integration method" << std::endl;
        throw std::logic_error("unkown integration method, abort program!");
    }
    return stream;
}

enum boundary_condition {
    no_boundaries,
    periodic
};

inline std::ostream &operator << (std::ostream &stream, const boundary_condition boundaries)
{
    if (boundaries == periodic)
        stream << "periodic boundaries" << std::endl;
    else if (boundaries == no_boundaries)
        stream << "no boundaries = open box" << std::endl;
    else {
        std::cerr << "unkown boundary condition" << std::endl;
        throw std::logic_error("unkown boundary condition, abort program!");
    }
    return stream;
}

enum lattice_types {
    onlyone,
    rectangular,
    triangular,
    random_lattice,
    individual
};

inline std::ostream &operator << (std::ostream &stream, const lattice_types lattice)
{
    if (lattice == onlyone)
        stream << "onlyone lattice" << std::endl;
    else if (lattice == rectangular)
        stream << "rectangular lattice" << std::endl;
    else if (lattice == triangular)
        stream << "triangular lattice" << std::endl;
    else if (lattice == random_lattice)
        stream << "random positioning of the particles" << std::endl;
    else if (lattice == individual)
        stream << "individual configuration chosen" << std::endl;
    else {
        std::cerr << "unkown lattice type" << std::endl;
        throw std::logic_error("unkown lattice type, abort program!");
    }
    return stream;
}
struct configuration{
    unsigned int number_particles;  ///number of particles in the simulation
    unsigned int seed;              ///seed of the random number generator
    double temperature;             ///temperature of the system
    double sigma;                   ///lengthscale of the potential
    double epsilon;                 ///strength of the potential
    double normed_distance;         ///initial distance between two neighboured particles (Gleichgewicht)
    double r_cut;                   ///cutting distance of the potential
    double dt;                      ///time discretization
    double box_width;               ///width  of the  simulation area
    double box_height;              ///height of the simulation area
    double space_in_x;              ///space between border and last particle
    double space_in_y;              ///space between border and last particle
    double axial_ratio;             ///ratio between the vertical and horizontal lattice rescale
    integration_method int_method;  ///integration method
    lattice_types lattice;          ///lattice, on which the particles are initialized. Choose between rectangular, triangular, random_lattice or individual configuration
    boundary_condition boundaries;  ///boundary_conditions, either no_boundaries or periodic
	double m_latticeConstant;       ///scaling of lattice
	double m_particleSpeed;         ///initial speed of particles
	bool    nogui;                  ///don't show graphical representation
	unsigned int    runs;           ///number of simulations steps
	std::string logname;            ///filename of the logfile
};

//overloaded << operator to show easily the configuration on the screen
inline std::ostream &operator << (std::ostream &stream, const configuration& my_config)
{
    stream << std::endl;
    stream << std::endl;
    stream << "chosen configuration:" << std::endl;
    stream << std::setw(21) << "number of particles: " << my_config.number_particles << std::endl;
    stream << std::endl;
    stream << std::setw(21) << "seed: " << my_config.seed << std::endl;
    stream << std::endl;
    stream << std::setw(21) << "temperature: " << my_config.temperature << std::endl;
    stream << std::endl;
    stream << std::setw(21) << "integration method: " << my_config.int_method << std::endl;
    stream << std::setw(21) << "timesteps dt: " << my_config.dt << std::endl;
    stream << std::endl;
    stream << "settings of the potential" << std::endl;
    stream << std::setw(21) << "lengthscale sigma: " << my_config.sigma << std::endl;
    stream << std::setw(21) << "normed distance between two neighboured particles:  " << my_config.normed_distance << std::endl;
    stream << std::setw(21) << "cutting distance of the potential (0.0 = infinity): " << my_config.r_cut << std::endl;
    stream << std::setw(21) << "strength epsilon of the potential: " << my_config.epsilon << std::endl;
    stream << std::endl;
    stream << "dimensions of the box (width x height): (" << my_config.box_width << " x " << my_config.box_height << ")" << std::endl;
    stream << "chosen boundary conditions: " << my_config.boundaries << std::endl;
    stream << std::endl;
    stream << std::setw(21) << "initial speed: " << my_config.m_particleSpeed << std::endl;
    stream << std::setw(21) << "lattice type:  " << my_config.lattice << std::endl;
    stream << std::setw(21) << "axial ratio:   " << my_config.axial_ratio << std::endl;
	stream << std::setw(21) << "lattice constant:   " << my_config.m_latticeConstant << std::endl;
	stream << std::setw(21) << "particle speed:   " << my_config.m_particleSpeed << std::endl;
    stream << std::endl;

    return stream;
}

#endif
