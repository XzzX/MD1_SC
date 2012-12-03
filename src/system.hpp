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
        configuration& m_config;

        ///list of kinetic energy of the whole system for every time step
        std::list<double> kin_energy;
        ///list of potential energy of the whole system for every time step
        std::list<double> pot_energy;
        ///trajectories of 10 randomly choosen particles used for diffusion calculation
        std::list<Vector> positionParticleList[10];
        ///position of the center of mass for every time step
        std::list<Vector> position_list;
        ///velocity of the center of mass for every time step
        std::list<Vector> velocity_list;
        //use these lists for storing position and velocity

        ///list of all simulated particles
        std::vector<particles> m_particles;
        ///saved starting positions of particles for diffusion calculation
        std::vector<Vector> m_startPosition;
        ///linear diffusion
        std::list<double> m_D;
        ///quadratic diffusion
        std::list<double> m_D2;
        ///temperature
        std::list<double> m_T;

        ///Lagrange-Multiplicator for constant temperature
        double  alpha;

        ///fricton constant for oscillator potential
        double  m_R;
        ///total time the system has evolved
        double  m_timeElapsed;

        ///handels subdivision of simulation box
		CellSubdivision*	m_CellSubdivision;

        ///returns number of particles which are simulated
        inline unsigned int GetNumberParticles() const {return m_particles.size();};

        ///constructor of the class
        molecular_dynamics(configuration &init_config);
        ///move the simulation one timestep dt
        void move_timestep();
        ///dump all saved observables to a file
        void    DumpData(const std::string& filename);

        ///seed the random generator with the given value
		void	RandomSeed(unsigned int seed);
		///generate a uniform distributed random number
		double	RandomUniform(double start, double stop);

        ///corrects the distance between two particles in correspondence to boundary conditions
		void	CorrectDistance(Vector& dist);
		void	CorrectPosition(particles& pos);

        ///save all observables
        void    Observe();

        ///calculates the correlation
        void    CalculateEndStats();

        ///updates the cell of the given particle
        void    UpdateParticleWithinCellSubdivision(const unsigned int i); //inline

		///returns the current temperature of the system
		double GetT();

    private:

        ///use this function to calculate the starting configuration of the particles
        void initialise_particles();
        ///initialise only one particle
        void InitParticlesOne();
        ///initialise particles in a rectangular shape
        void InitParticlesRectangular();
        ///initalise particles in a triangular shape
        void InitParticlesTriangular();

		double potentialO(const Vector& pos, const double eps=0.5);
		double potentialLJ(const Vector& rij);
		Vector forceO(const Vector& pos, const double eps=0.5);
		Vector forceLJ(const Vector& rij);

        void IntegrateEuler();
        void IntegrateLeapFrog();

        ///use this function to calculate the acceleration of the particle
        void CalculateAcceleration();

        void SetCMSP0();

        ///updates lagrange mulitplicator alpha
        void UpdateAlpha();

        ///calculates the correlation
        void    Correlate();
        ///calculates the speed distribution
        void    VelocityDistribution();
};

inline
void    molecular_dynamics::UpdateParticleWithinCellSubdivision(const unsigned int i){
    m_CellSubdivision->DeleteParticle(m_particles[i].m_cell_id, i);
    m_particles[i].m_cell_id = m_CellSubdivision->InsertParticle(m_particles[i].m_position, i);
}

#endif
