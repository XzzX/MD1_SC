#include <iostream>
#include <list>
#include <fstream>
#include    <string>

#include "system.hpp"

#include    "Histogram.h"

using namespace std;

molecular_dynamics::molecular_dynamics(configuration &init_config) : m_generator(init_config.seed), m_R(0.0), m_timeElapsed(0.0), m_config(init_config){

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

/**
Calculates the minimum distance between two particles in correspondence to boundary conditions.
@param dist plain distance between two particles
**/
void	molecular_dynamics::CorrectDistance(Vector& dist){
	if (dist[0]>(m_config.box_width*0.5)) dist[0]-=m_config.box_width; else
		if (dist[0]<(-m_config.box_width*0.5)) dist[0]+=m_config.box_width;
	if (dist[1]>(m_config.box_height*0.5)) dist[1]-=m_config.box_height; else
		if (dist[1]<(-m_config.box_height*0.5)) dist[1]+=m_config.box_height;
}

/**
Calculates the position correction with paying attention to boundary conditions. Also saves the count of border crossings.
@param part particle which should be corrected
**/
void	molecular_dynamics::CorrectPosition(particles& part){
	if (part.m_position[0]>(m_config.box_width)){
	    part.m_position[0] -= m_config.box_width;
	    part.m_borderCrossingX++;
    }else if (part.m_position[0]<(0.0)){
        part.m_position[0] += m_config.box_width;
        part.m_borderCrossingX--;
    }
	if (part.m_position[1]>(m_config.box_height)){
	    part.m_position[1] -= m_config.box_height;
	    part.m_borderCrossingY++;
    }else if (part.m_position[1]<(0.0)){
        part.m_position[1] += m_config.box_height;
        part.m_borderCrossingY--;
    }
}

/**
Dumps observed data to file. Observed informations are:
time
position
speed
kinetic energy
potential energy
distance to start position
distance to start position ^2
@param filename filename of the dump file
**/
void    molecular_dynamics::DumpData(const std::string& filename){
    fstream fout(filename.c_str(), fstream::out);
    fstream fout2((std::string("rawpos_")).append(filename).c_str(), fstream::out);

	fout << m_config;
	fout2 << m_config;

    list<Vector>::iterator pos0 = positionParticleList[0].begin();
    list<Vector>::iterator pos1 = positionParticleList[1].begin();
    list<Vector>::iterator pos2 = positionParticleList[2].begin();
    list<Vector>::iterator pos3 = positionParticleList[3].begin();
    list<Vector>::iterator pos4 = positionParticleList[4].begin();
    list<Vector>::iterator pos5 = positionParticleList[5].begin();
    list<Vector>::iterator pos6 = positionParticleList[6].begin();
    list<Vector>::iterator pos7 = positionParticleList[7].begin();
    list<Vector>::iterator pos8 = positionParticleList[8].begin();
    list<Vector>::iterator pos9 = positionParticleList[9].begin();

    list<Vector>::iterator it0 = position_list.begin();
    list<Vector>::iterator it1 = velocity_list.begin();
    list<double>::iterator it2 = kin_energy.begin();
    list<double>::iterator it3 = pot_energy.begin();
    list<double>::iterator it4 = m_D.begin();
    list<double>::iterator it5 = m_D2.begin();
    for(unsigned int i=0; i<position_list.size(); i++, it0++, it1++, it2++, it3++, it4++, it5++, pos0++, pos1++, pos2++, pos3++, pos4++, pos5++, pos6++, pos7++, pos8++, pos9++){
        fout << i*m_config.dt << "\t";
        fout << *it0 << "\t";
        fout << *it1 << "\t";
		fout << *it2 << "\t";
		fout << *it3 << "\t";
		fout << *it4 << "\t";
        fout << *it5 << endl;

        fout2 << *pos0 << "\t" << *pos1 << "\t" << *pos2 << "\t" << *pos3 << "\t" << *pos4 << "\t" << *pos5 << "\t" << *pos6 << "\t" << *pos7 << "\t" << *pos8 << "\t" << *pos9 << "\t" << endl;
    }
    fout.close();
    fout2.close();
}

void molecular_dynamics::initialise_particles(){
    m_particles.clear();
    if (m_config.lattice == onlyone) InitParticlesOne();
    if (m_config.lattice == rectangular) InitParticlesRectangular();
    if (m_config.lattice == triangular) InitParticlesTriangular();

    SetCMSP0();

    m_config.number_particles = GetNumberParticles();
    m_startPosition.clear();
    for (unsigned int i=0; i<GetNumberParticles(); i++)
        m_startPosition.push_back(m_particles[i].m_position);
}

/**
Initialises one particle
**/
void molecular_dynamics::InitParticlesOne(){
    //initialise the particle
    particles dummy;

    EmptyParticle(dummy);
    dummy.m_position = Vector(20.0, 0.0, 0.0);
    dummy.m_speed = Vector(0.0, 5.0, 0.0);

    m_particles.push_back(dummy);
}

/**
Tries to initialise m_config.number_particles particles in a rectangular pattern.
**/
void molecular_dynamics::InitParticlesRectangular(){
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

            EmptyParticle(dummy);
            dummy.m_position = i*a + j*b + delta;
            double alpha = RandomUniform(0.0, 2*M_PI);
            dummy.m_speed = Vector(cos(alpha), sin(alpha), 0.0) * m_config.m_particleSpeed;
            dummy.m_radius = 0.2;

            m_particles.push_back(dummy);
        }

    m_config.box_height = (numb-1) * m_config.m_latticeConstant * m_config.sigma + 2*m_config.space_in_y;
    m_config.box_width  = (numa-1) * m_config.m_latticeConstant * m_config.sigma + 2*m_config.space_in_x;
}

/**
Tries to initialise m_config.number_particles particles in a triangular pattern.
**/
void molecular_dynamics::InitParticlesTriangular(){
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

            EmptyParticle(dummy);
            dummy.m_position = i*a + j*b - (j/2)*a + delta;
            double alpha = RandomUniform(0.0, 2*M_PI);
            dummy.m_speed = Vector(cos(alpha), sin(alpha), 0.0) * m_config.m_particleSpeed;
            dummy.m_radius = 0.2;

            m_particles.push_back(dummy);
        }

    m_config.box_height = b[1] * (numb-1) + 2*m_config.space_in_y;
    m_config.box_width  = (numa-1) * m_config.m_latticeConstant * m_config.sigma + 2*m_config.space_in_x;
}

/**
Calculates the potential of the harmonic oscillator 0.5*eps*x^2
@param pos position of the particle
@param eps parameter of the harmonic oscillator potential
@return value of the potential
**/
double molecular_dynamics::potentialO(const Vector& pos, const double eps){
    return -0.5*eps*norm(pos);
}

/**
Calculates the Lennard-Jones-Potential
@param rij distance from the particle
@return value of the potential
**/
double molecular_dynamics::potentialLJ(const Vector& rij){
	//Lenard-Jones
	double r = norm(rij);
	if (r>m_config.r_cut)
        return 0.0;
	else
        return 4.0 * m_config.epsilon * ( pow(m_config.sigma/r, 12) - pow(m_config.sigma/r,6) ) + m_config.epsilon;
}

/**
Calculates the resulting force of potential of the harmonic oscillator 0.5*eps*x^2
@param pos position of the particle
@param eps parameter of the harmonic oscillator potential
@return value of the potential
**/
Vector molecular_dynamics::forceO(const Vector& pos, const double eps){
	//harmonic oscillator
    return -eps*pos;
}

/**
Calculates the resulting force of Lennard-Jones-Potential
@param rij distance from the particle
@return value of the potential
**/
Vector molecular_dynamics::forceLJ(const Vector& rij){
	//Lenard-Jones
	double r = norm(rij);
	if (r>m_config.r_cut)
        return Vector(0.0,0.0,0.0);
	else
        return 24.0 * m_config.epsilon / m_config.sigma / m_config.sigma * ( pow(m_config.sigma/r, 8) - 2 * pow(m_config.sigma/r,14) ) * rij;
}

void molecular_dynamics::move_timestep()
{
    if (m_config.int_method==euler) IntegrateEuler();
    else if (m_config.int_method==leap_frog) IntegrateLeapFrog();

    m_timeElapsed += m_config.dt;
}

/**
Integrates the equation of motion with Euler's method
**/
void molecular_dynamics::IntegrateEuler(){
    for(unsigned int i=0; i<GetNumberParticles(); i++){
        m_particles[i].m_position += m_particles[i].m_speed*m_config.dt;
        m_particles[i].m_speed += m_particles[i].m_acceleration*m_config.dt*1.0;
        CalculateAcceleration();
		if (m_config.boundaries==periodic){
		    CorrectPosition(m_particles[i]);
            UpdateParticleWithinCellSubdivision(i);
		}
		//m_particles[i].m_color = m_CellSubdivision->GetCellColor(m_particles[i].m_cell_id);
	}
}

/**
Integrates the equation of motion with LeapFrog method
**/
void molecular_dynamics::IntegrateLeapFrog(){
    for(unsigned int i=0; i<GetNumberParticles(); i++){

		m_particles[i].m_speed += m_particles[i].m_acceleration*m_config.dt*0.5;
        m_particles[i].m_position += m_particles[i].m_speed*m_config.dt;
	}
	CalculateAcceleration();
	for(unsigned int i=0; i<GetNumberParticles(); i++){
		m_particles[i].m_speed += m_particles[i].m_acceleration*m_config.dt*0.5;

		if (m_config.boundaries==periodic){
		    CorrectPosition(m_particles[i]);
            UpdateParticleWithinCellSubdivision(i);
		}
		//m_particles[i].m_color = m_CellSubdivision->GetCellColor(m_particles[i].m_cell_id);
	}
}

/**
Calculates the actual acceleration within the current potential
**/
void molecular_dynamics::CalculateAcceleration(){
    list<unsigned int>::iterator    it;
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

void molecular_dynamics::UpdateAlpha(){
    alpha = 0.0;
    double alpha2(0.0);

    for (unsigned int i=0; i<GetNumberParticles(); i++){
        alpha -= dot(m_particles[i].m_position, m_particles[i].m_acceleration);
        alpha2 += dot(m_particles[i].m_speed, m_particles[i].m_speed);
    }

    alpha /= alpha2;
}

void molecular_dynamics::Observe(){
    if ((m_timeElapsed/m_config.dt)<100000){
        Vector pos(0,0,0);
        Vector v(0,0,0);
        double  Ekin(0);
        double  Epot(0);
        double  D = 0.0;
        double  D2 = 0.0;
        Vector  a(m_config.box_width,0.0,0.0);
        Vector  b(0.0,m_config.box_height,0.0);
        for (unsigned int i = 0; i < GetNumberParticles(); i++){
            pos += m_particles[i].m_position;
			v += m_particles[i].m_speed;
			Ekin += 0.5*m_particles[i].m_mass*norm_sq(m_particles[i].m_speed);
			for (unsigned int j = i+1; j < GetNumberParticles(); j++){
			    Vector rij(m_particles[j].m_position - m_particles[i].m_position);
			    if (m_config.boundaries==periodic) CorrectDistance(rij);
				Epot += potentialLJ(rij);
			}
			double  d = norm(m_particles[i].m_position+m_particles[i].m_borderCrossingX*a+m_particles[i].m_borderCrossingY*b - m_startPosition[i]);
			D += d;
			D2 += d*d;
        }
        position_list.push_back(pos / GetNumberParticles());
        velocity_list.push_back(v / GetNumberParticles());
        pot_energy.push_back(Epot);
        kin_energy.push_back(Ekin);
        m_D.push_back(D/GetNumberParticles());
        m_D2.push_back(D2/GetNumberParticles());

        positionParticleList[0].push_back(m_particles[111].m_position);
        positionParticleList[1].push_back(m_particles[222].m_position);
        positionParticleList[2].push_back(m_particles[333].m_position);
        positionParticleList[3].push_back(m_particles[444].m_position);
        positionParticleList[4].push_back(m_particles[555].m_position);
        positionParticleList[5].push_back(m_particles[666].m_position);
        positionParticleList[6].push_back(m_particles[777].m_position);
        positionParticleList[7].push_back(m_particles[888].m_position);
        positionParticleList[8].push_back(m_particles[561].m_position);
        positionParticleList[9].push_back(m_particles[134].m_position);
    }
}

void    molecular_dynamics::Correlate(){
    //histogram
    static const unsigned int bins = 500;

    Histogram   histo(0.0, m_config.box_width*0.5, bins);

    fstream fout(std::string("corraw_").append(m_config.logname).c_str(), fstream::out);

	fout << m_config;

    for (unsigned int i = 0; i < GetNumberParticles(); i++)
        for (unsigned int j = i+1; j < GetNumberParticles(); j++){
            Vector rij(m_particles[j].m_position - m_particles[i].m_position);
            if (m_config.boundaries==periodic) CorrectDistance(rij);
            histo.Add(norm(rij));
            fout << norm(rij) << std::endl;
        }

    fout.close();


    //dump
    fout.open(std::string("cor_").append(m_config.logname).c_str(), fstream::out);
	fout << m_config;
    for (unsigned int i = 0; i < bins; i++)
        fout << histo.GetXUpper(i) <<  "\t" << histo.GetValue(i) <<  "\t" <<
        (m_config.box_height*m_config.box_width*histo.GetValue(i)/2.0/M_PI/histo.GetBinWidth()/histo.GetBinWidth()/m_config.number_particles/m_config.number_particles/(static_cast<double>(i)-0.5)) << std::endl;

    fout.close();

    std::cout << "box_height: " << m_config.box_height << std::endl;
    std::cout << "box_width : " << m_config.box_height << std::endl;
    std::cout << "N         : " << m_config.number_particles << std::endl;
}

void molecular_dynamics::SetCMSP0(){
    //berechne Gesamtimpuls
    Vector p(0,0,0);

    for (unsigned int i=0; i<m_particles.size(); i++){
        p += m_particles[i].m_speed * m_particles[i].m_mass;
    }

    p /= static_cast<double>(m_particles.size());

    std::cout << "Gesamtimpuls: " << p << std::endl;

    //Neuberechnung der Impulse

    for (unsigned int i=0; i<m_particles.size(); i++){
        m_particles[i].m_speed -= p / m_particles[i].m_mass;
    }
}

void molecular_dynamics::CalculateEndStats(){
    Correlate();
    VelocityDistribution();
}

void molecular_dynamics::VelocityDistribution(){
    fstream fout((std::string("speed_")).append(m_config.logname).c_str(), fstream::out);
	fout << m_config;

	for (unsigned int i=0; i<GetNumberParticles(); i++){
        fout << norm(m_particles[i].m_speed) << std::endl;
	}

	fout.close();
}
