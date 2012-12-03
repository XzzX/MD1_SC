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
        stream << "eulerstep";
    else if (method == leap_frog)
        stream << "leapfrog integration";
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
        stream << "periodic boundaries";
    else if (boundaries == no_boundaries)
        stream << "no boundaries = open box";
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
        stream << "onlyone lattice";
    else if (lattice == rectangular)
        stream << "rectangular lattice";
    else if (lattice == triangular)
        stream << "triangular lattice";
    else if (lattice == random_lattice)
        stream << "random positioning of the particles";
    else if (lattice == individual)
        stream << "individual configuration chosen";
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
	stream << "#CONFIGURATION" << std::endl;
	stream << "#number of particles: " << my_config.number_particles << std::endl;
	stream << "#seed:                " << my_config.seed << std::endl;
	stream << "#temperature:         " << my_config.temperature << std::endl;
	stream << "#normed distance:     " << my_config.normed_distance << std::endl;
	stream << "#r_cut:               " << my_config.r_cut << std::endl;
	stream << "#dt:                  " << my_config.dt << std::endl;
	stream << "#box_width:           " << my_config.box_width << std::endl;
	stream << "#box height:          " << my_config.box_height << std::endl;
	stream << "#space_in_x:          " << my_config.space_in_x << std::endl;
	stream << "#space_in_y:          " << my_config.space_in_y << std::endl;
	stream << "#axial_ratio:         " << my_config.axial_ratio << std::endl;
	stream << "#integration method:  " << my_config.int_method << std::endl;
	stream << "#lattice type:        " << my_config.lattice << std::endl;
	stream << "#boundary conditions: " << my_config.boundaries << std::endl;
	stream << "#lattice constant:    " << my_config.m_latticeConstant << std::endl;
	stream << "#particle speed:      " << my_config.m_particleSpeed << std::endl;
	stream << "#nogui:               " << my_config.nogui << std::endl;
	stream << "#runs:                " << my_config.runs << std::endl;
	stream << "#logname:             " << my_config.logname << std::endl;
    return stream;
}

#endif
