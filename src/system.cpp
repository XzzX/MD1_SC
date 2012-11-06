#include <iostream>
#include <list>
#include <fstream>

#include "system.hpp"

#include    "Histogram.h"

using namespace std;

molecular_dynamics::molecular_dynamics(const configuration &init_config) : m_timeElapsed(0.0), m_R(0.0), m_generator(init_config.seed){
	(this)->m_config = init_config;

    initialise_particles();

	m_CellSubdivision = new CellSubdivision(m_config.box_width, m_config.box_height, m_config.r_cut);
	for (unsigned int i=0; i<GetNumberParticles(); i++){
		m_particles[i].m_cell_id = m_CellSubdivision->InsertParticle(m_particles[i].m_position, i);
		m_particles[i].m_color = m_CellSubdivision->GetCellColor(m_particles[i].m_cell_id);
	}
}

// http://www.boost.org/doc/libs/1_51_0/libs/random/example/random_demo.cpp

void	molecular_dynamics::RandomSeed(unsigned int seed){
	m_generator.seed(seed);
}
double	molecular_dynamics::RandomUniform(double start, double stop){
	boost::uniform_real<> uni_dist(start,stop);
	boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(m_generator, uni_dist);
	return uni();
}

void	molecular_dynamics::CorrectDistance(Vector& dist){
	if (dist[0]>(m_config.box_width*0.5)) dist[0]-=m_config.box_width; else
		if (dist[0]<(-m_config.box_width*0.5)) dist[0]+=m_config.box_width;
	if (dist[1]>(m_config.box_height*0.5)) dist[1]-=m_config.box_height; else
		if (dist[1]<(-m_config.box_height*0.5)) dist[1]+=m_config.box_height;
}
void	molecular_dynamics::CorrectPosition(Vector& pos){
	if (pos[0]>(m_config.box_width))	pos[0]-=m_config.box_width; else
		if (pos[0]<(0.0))					pos[0]+=m_config.box_width;
	if (pos[1]>(m_config.box_height))	pos[1]-=m_config.box_height; else
		if (pos[1]<(0.0))					pos[1]+=m_config.box_height;
}

void    molecular_dynamics::DumpData(const std::string& filename){
    fstream fout(filename.c_str(), fstream::out);

    list<Vector>::iterator it0 = position_list.begin();
    list<Vector>::iterator it1 = velocity_list.begin();
    list<double>::iterator it2 = kin_energy.begin();
    list<double>::iterator it3 = pot_energy.begin();
    for(unsigned int i=0; i<position_list.size(); i++, it0++, it1++, it2++, it3++){
        fout << i*m_config.dt << "\t";
        fout << *it0 << "\t";
        fout << *it1 << "\t";
		fout << *it2 << "\t";
        fout << *it3 << endl;

    }
    fout.close();
}

void molecular_dynamics::initialise_particles(){
	std::cout << "initialising particles" << std::endl;
    m_particles.clear();
    if (m_config.lattice == onlyone){
        //initialise the particle
        particles dummy;
        dummy.m_position = Vector(20.0, 0.0, 0.0);
        dummy.m_speed = Vector(0.0, 5.0, 0.0);
        dummy.m_acceleration = Vector(0.0, 0.0, 0.0);
        dummy.m_radius = 2.0;
        dummy.m_mass = 1.0;
        dummy.m_color = c_white;

        m_particles.push_back(dummy);
		std::cout << "1 particle initialised" << std::endl;
		m_config.number_particles = 1;
    } else
    if (m_config.lattice == rectangular){
        //axial_ratio = v/h
		Vector	a(1.0, 0, 0);
		Vector	b(0, 1.0, 0);

		a *= m_config.m_latticeConstant * m_config.sigma;
		b *= m_config.m_latticeConstant * m_config.sigma;

		a *= m_config.sigma;
		b *= m_config.sigma;

		unsigned int	numa(floor(sqrt(m_config.number_particles/m_config.axial_ratio)));
		unsigned int	numb(floor(numa * m_config.axial_ratio));

		Vector delta(m_config.space_in_x,m_config.space_in_y, 0.0);

		for (unsigned int i = 0; i<numa; i++)
			for (unsigned int j = 0; j<numb; j++){
				particles dummy;
				dummy.m_position = i*a + j*b + delta;
				double alpha = RandomUniform(0.0, 2*M_PI);
				dummy.m_speed = Vector(cos(alpha), sin(alpha), 0.0) * m_config.m_particleSpeed;
				dummy.m_acceleration = Vector(0.0, 0.0, 0.0);
				dummy.m_radius = 0.2;
				dummy.m_mass = 1.0;
				dummy.m_color = c_white;

				m_particles.push_back(dummy);
			}

		m_config.box_height = (numb-1) * m_config.m_latticeConstant * m_config.sigma + 2*m_config.space_in_y;
		m_config.box_width  = (numa-1) * m_config.m_latticeConstant * m_config.sigma + 2*m_config.space_in_x;

		std::cout << numa*numb <<" particles initialised" << std::endl;
		m_config.number_particles = numa*numb;
    } else
    if (m_config.lattice == triangular){
        //axial_ratio = v/h
		Vector	a(1.0, 0, 0);
		Vector	b(cos(M_PI/3.0), sin(M_PI/3.0), 0);

		a *= m_config.m_latticeConstant * m_config.sigma;
		b *= m_config.m_latticeConstant * m_config.sigma;

		unsigned int	numa(floor(sqrt(m_config.number_particles/m_config.axial_ratio)));
		unsigned int	numb(floor(numa * m_config.axial_ratio));

		Vector delta(m_config.space_in_x,m_config.space_in_y, 0.0);

		for (unsigned int i = 0; i<numa; i++)
			for (unsigned int j = 0; j<numb; j++){
				particles dummy;
				dummy.m_position = i*a + j*b - (j/2)*a + delta;
				double alpha = RandomUniform(0.0, 2*M_PI);
				dummy.m_speed = Vector(cos(alpha), sin(alpha), 0.0) * m_config.m_particleSpeed;
				dummy.m_acceleration = Vector(0.0, 0.0, 0.0);
				dummy.m_radius = 0.2;
				dummy.m_mass = 1.0;
				dummy.m_color = c_white;

				m_particles.push_back(dummy);
			}

		m_config.box_height = b[1] * (numb-1) + 2*m_config.space_in_y;
		m_config.box_width  = (numa-1) * m_config.m_latticeConstant * m_config.sigma + 2*m_config.space_in_x;

		std::cout << numa*numb <<" particles initialised" << std::endl;
		m_config.number_particles = numa*numb;
    }
}

double molecular_dynamics::potentialO(const Vector& pos, const double eps){
    return -0.5*eps*norm(pos);
}

double molecular_dynamics::potentialLJ(const Vector& rij){
    return -0.5*norm(rij);
}

Vector molecular_dynamics::forceO(const Vector& pos, const double eps){
	//harmonic oscillator
    return -eps*pos;
}

Vector molecular_dynamics::forceLJ(const Vector& rij){
	//Lenard-Jones
	double r = norm(rij);
	return 24.0 * m_config.epsilon / m_config.sigma / m_config.sigma * ( pow(m_config.sigma/r, 8) - 2 * pow(m_config.sigma/r,14) ) * rij;
}

void molecular_dynamics::move_timestep()
{
    integrate_equation_of_motion();
    m_timeElapsed += m_config.dt;
}

void molecular_dynamics::integrate_equation_of_motion(){
    for(unsigned int i=0; i<GetNumberParticles(); i++){

		m_particles[i].m_speed += m_particles[i].m_acceleration*m_config.dt*0.5;
        m_particles[i].m_position += m_particles[i].m_speed*m_config.dt;
		//m_particles[i].m_speed += m_particles[i].m_acceleration*m_config.dt*1.0;
	}
	calculate_acceleration();
	for(unsigned int i=0; i<GetNumberParticles(); i++){
		m_particles[i].m_speed += m_particles[i].m_acceleration*m_config.dt*0.5;

		if (m_config.boundaries==periodic) CorrectPosition(m_particles[i].m_position);
		m_CellSubdivision->DeleteParticle(m_particles[i].m_cell_id, i);
		m_particles[i].m_cell_id = m_CellSubdivision->InsertParticle(m_particles[i].m_position, i);
		//m_particles[i].m_color = m_CellSubdivision->GetCellColor(m_particles[i].m_cell_id);
	}
}

void molecular_dynamics::calculate_acceleration(){
    list<int>::iterator    it;
	for(unsigned int i=0; i<GetNumberParticles(); i++){

		Vector F(0,0,0);
		//F += force(m_particles[i].m_position);
		//F += m_R * m_particles[i].m_speed;
		/*for(unsigned int j=0; j<GetNumberParticles(); j++)
			if (i!=j) {
				Vector rij(m_particles[j].m_position - m_particles[i].m_position);
				if (m_config.boundaries==periodic) CorrectDistance(rij);

				F += forceLJ(rij);
			}*/

        m_CellSubdivision->GetNeighbours(m_particles[i].m_cell_id, i);
        for (it=m_CellSubdivision->GetNeighbourBegin(); it!=m_CellSubdivision->GetNeighbourEnd(); it++){
            Vector rij(m_particles[*it].m_position - m_particles[i].m_position);
            if (m_config.boundaries==periodic) CorrectDistance(rij);
            F += forceLJ(rij);
        }
        /*m_particles[i].m_color= c_white;
        if (i==960){
            m_particles[i].m_color= c_green;
            for (it=m_CellSubdivision->GetNeighbourBegin(); it!=m_CellSubdivision->GetNeighbourEnd(); it++)
                m_particles[*it].m_color = c_red;
        }*/

		m_particles[i].m_acceleration = F/m_particles[i].m_mass;
	}
}

void molecular_dynamics::Observe(){
    if ((m_timeElapsed/m_config.dt)<100000){
        Vector pos(0,0,0);
        Vector v(0,0,0);
        double  Ekin(0);
        double  Epot(0);
        for (unsigned int i = 0; i < GetNumberParticles(); i++){
            pos += m_particles[i].m_position;
			v += m_particles[i].m_speed;
			Ekin += 0.5*m_particles[i].m_mass*norm_sq(m_particles[i].m_speed);
			for (unsigned int j = i+1; j < GetNumberParticles(); j++){
			    Vector rij(m_particles[j].m_position - m_particles[i].m_position);
			    if (m_config.boundaries==periodic) CorrectDistance(rij);
				Epot += potentialLJ(rij);
			}
        }
        position_list.push_back(pos / GetNumberParticles());
        velocity_list.push_back(v / GetNumberParticles());
        pot_energy.push_back(Epot);
        kin_energy.push_back(Ekin);
    }
}

void    molecular_dynamics::Correlate(){
    //histogram
    static const unsigned int bins = 200;

    Histogram   histo(0.0, m_config.box_width*0.5, bins);

    for (unsigned int i = 0; i < GetNumberParticles(); i++)
        for (unsigned int j = i+1; j < GetNumberParticles(); j++){
            Vector rij(m_particles[j].m_position - m_particles[i].m_position);
            if (m_config.boundaries==periodic) CorrectDistance(rij);
            histo.Add(norm(rij));
        }


    //dump
    fstream fout("correlation.txt", fstream::out);

    for (unsigned int i = 0; i < bins; i++)
        fout << histo.GetXUpper(i) <<  "\t" << histo.GetValue(i) <<  "\t" <<
        (m_config.box_height*m_config.box_width*histo.GetValue(i)/2.0/M_PI/histo.GetBinWidth()/histo.GetBinWidth()/m_config.number_particles/m_config.number_particles/(static_cast<double>(i)-0.5)) << std::endl;

    fout.close();
}
