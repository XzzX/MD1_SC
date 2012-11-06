#ifndef SYSTEM_HPP
#define SYSTEM_HPP
#include <sstream>
//#include <omp.h>

#include <cstring>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <list>
#include <limits>

#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

#include "Vector.hpp"
#include "Window.hpp"
#include "configuration.hpp"
#include "particles.hpp"
#include "CellSubdivision.h"

typedef boost::minstd_rand base_generator_type;

class molecular_dynamics {
    public:
		///random number generator
		base_generator_type		m_generator;

        ///configuration structure with all parameters of the system
        configuration m_config;

        std::list<double> kin_energy;
        std::list<double> pot_energy;
        //use these lists to store the values of kinetic and potential energy

        std::list<Vector> position_list;
        std::list<Vector> velocity_list;
        //use these lists for storing position and velocity

        std::vector<particles> m_particles;

        double  m_R;
        double  m_timeElapsed;

		CellSubdivision*	m_CellSubdivision;

        inline unsigned int GetNumberParticles() const {return m_particles.size();};

        ///constructor of the class
        molecular_dynamics(const configuration &init_config);
        ///move the simulation one timestep dt
        void move_timestep();
        ///dump all saved observables to a file
        void    DumpData(const std::string& filename);

		void	RandomSeed(unsigned int seed);
		double	RandomUniform(double start, double stop);

		void	CorrectDistance(Vector& dist);
		void	CorrectPosition(Vector& pos);

        ///save all observables
        void    Observe();

        void    Correlate();

    private:

        ///use this function to calculate the starting configuration of the particles
        void initialise_particles();

		double potentialO(const Vector& pos, const double eps=0.5);
		double potentialLJ(const Vector& rij);
		Vector forceO(const Vector& pos, const double eps=0.5);
		Vector forceLJ(const Vector& rij);

        ///this function integrates the equation of motion either by using euler steps or by leapfrog integration
        void integrate_equation_of_motion();

        ///use this function to calculate the acceleration of the particle
        void calculate_acceleration();
};

#endif
