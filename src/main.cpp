#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <boost/timer/timer.hpp>

#include "Window.hpp"
#include "Vector.hpp"
#include "system.hpp"
#include "conversions.hpp"
#include "configuration.hpp"
#include "particles.hpp"
#include "plot_utils.hpp"
#include "MeanVar.h"

//installation unter Ubuntu 10.04
//sudo apt-get install g++
//sudo apt-get install libboost-dev libboost-doc
//download sfml unter http://www.sfml-dev.org/download.php
//entpacke den sfml Ordner
//Im Ordner sfml-x.y sudo make install eingeben
//sudo apt-get install freeglut3 freeglut3-dev libglew1.5 libglew1.5-dev

int main( int argc, char **argv ) {

    configuration my_config;

    my_config.seed = 0;
    my_config.number_particles = 0;
    my_config.int_method = leap_frog;



    if (argc == 1) //standard configuration
    {
        std::cout << "standard configuration" << std::endl;
        my_config.seed = 49;
        my_config.number_particles = 1000;
        my_config.sigma = 1.0;
        my_config.epsilon = 1.0;
        my_config.temperature = 1.0;
        my_config.normed_distance = 1.0*pow(2.0,1.0/6.0); //(Gleichgewicht)
//        my_config.r_cut  = 2.5 * my_config.sigma;
        my_config.r_cut  = pow(2.0,1.0/6.0) * my_config.sigma;
        my_config.int_method = leap_frog;
        my_config.dt = 0.01;
        my_config.v_initial = 0.;
        my_config.lattice = triangular;
        my_config.axial_ratio = 1.0; // v/h
        my_config.box_width = 100.0;
        my_config.box_height = 100.0;
        my_config.boundaries = periodic;
		my_config.m_latticeConstant = my_config.normed_distance * 1.1;
        my_config.space_in_x = my_config.m_latticeConstant * 0.5;
        my_config.space_in_y = my_config.m_latticeConstant * 0.5;
		my_config.m_particleSpeed = 1.0;
    } else { //individual configuration
        for (int i = 1; i < argc-1; i += 2) //parse input parameters
        {
            std::string my_string = argv[i];
            if (my_string.compare("-N") == 0)
            {
                my_config.number_particles = string_to_int(argv[i+1]);
            }
            else if (my_string.compare("-seed") == 0)
            {
                my_config.seed = string_to_int(argv[i+1]);
            }
            else if (my_string.compare("-int_method") == 0)
            {
                my_config.int_method = string_to_integration_method(argv[i+1]);
            }
        }
    }

    const Vector main_size(800,600);
    const Vector secondary_size(400,300);

    Window main_frame("Lennard Jones Fluid", main_size, 10);

	std::cout << my_config << std::endl;

    molecular_dynamics md_system(my_config);

	main_frame.SetCameraPosition(Vector(-md_system.m_config.box_width*0.5, -md_system.m_config.box_height*0.5, -md_system.m_config.box_height*0.5));

	MeanVar	perf;

	boost::timer::cpu_timer timer;

    while (main_frame.isOpen()){
    	main_frame.clear();
    	main_frame.processEvents();
		md_system.Observe();
		timer.start();
        md_system.move_timestep();
        timer.stop();
        perf.Add(timer.elapsed().user);
        if (md_system.m_config.boundaries==periodic){
			Vector pos(md_system.m_config.box_width*0.5, md_system.m_config.box_height*0.5, 0.0);
			main_frame.drawRectangle(pos,md_system.m_config.box_width,md_system.m_config.box_height,c_white);
		}
        for(unsigned int i=0; i<md_system.GetNumberParticles(); i++)
            drawParticle(md_system.m_particles[i], main_frame);
 	    main_frame.display();
    }
	std::cout << perf.Mean() << std::endl;
	md_system.Correlate();
    md_system.DumpData("data.txt");
    return 0;
}
