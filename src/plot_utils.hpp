#ifndef PLOT_UTILS_HPP
#define PLOT_UTILS_HPP
#include "Window.hpp"
#include "particles.hpp"
#include <list>
#include <vector>
#include <limits>

struct arrangement {
    double x_range_min;
    double x_range_max;

    double y_range_min;
    double y_range_max;

    bool fixed;
};

void drawParticle(const particles &my_particle, const Window &w);
void drawCellSubdivision(const unsigned int n_horizontal, const unsigned int n_vertical, const double width, const double height, const cColor& color, const Window &w);
void drawValues(const std::pair<std::list<double>::const_iterator, std::list<double>::const_iterator> &values, arrangement &pos, const cColor& color, const Window &w);
void drawCoordinateAxis(arrangement &pos, const cColor& color, const Window &w);


inline void drawParticle(const particles &my_particle, const Window &w)
{
    w.drawCircle(my_particle.m_position, my_particle.m_color, my_particle.m_radius);
}


#endif
