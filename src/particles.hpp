#ifndef PARTICLES_HPP
#define PARTICLES_HPP
#include "Vector.hpp"
#include <iomanip>

struct cColor{
    double r;
    double g;
    double b;
};

const cColor c_red   = {1.0f, 0.0f, 0.0f};
const cColor c_green = {0.0f, 1.0f, 0.0f};
const cColor c_blue  = {0.0f, 0.0f, 1.0f};
const cColor c_white = {1.0f, 1.0f, 1.0f};

inline bool periodic_vector(Vector &input_vector, const double width, const double height)
{
    //periodic vector returns per reference a new vector with respect to periodic boundaries of a given box
    //with dimensions width and height. If the vector is changed, it would return true
    bool is_periodic = false;
    if (input_vector[0] > width / 2.0)
    {
        while (input_vector[0] > width / 2.0)
            input_vector[0] -= width;

        is_periodic = true;
    } else if (input_vector[0] <= -width / 2.0)
    {
        while (input_vector[0] <= -width / 2.0)
            input_vector[0] += width;
        is_periodic = true;
    }

    if (input_vector[1] > height / 2.0)
    {
        while (input_vector[1] > height / 2.0)
            input_vector[1] -= height;
        is_periodic = true;
    } else if (input_vector[1] <= -height / 2.0)
    {
        while (input_vector[1] <= -height / 2.0)
            input_vector[1] += height;
        is_periodic = true;
    }

    return is_periodic;
}

inline std::ostream &operator << (std::ostream &stream, const cColor& color)
{
    //overloaded ostream << operator to display the rgb values of a color
    stream << "(" << color.r << " , " << color.g << " , " << color.b << std::endl;
    return stream;
}

cColor createColor(unsigned int i, unsigned int N);

inline void compoment_to_color(double r, double g, double b, cColor &color)
{
    color.r = r;
    color.g = g;
    color.b = b;
}

struct particles {
    //structure which contains all the necessary informations about the particle
    Vector m_position;
    Vector m_speed;
    Vector m_acceleration;
    double m_radius;
    double m_mass;
    cColor m_color;
    unsigned int m_cell_id;
};

inline std::ostream &operator << (std::ostream &stream, const particles& particle)
{
    //overloaded ostream << operator to display the current information about one single particle
    stream << std::setw(20) << "particle position: "     << particle.m_position << std::endl;
    stream << std::setw(20) << "speed: "        << particle.m_speed << std::endl;
    stream << std::setw(20) << "acceleration: " << particle.m_acceleration << std::endl;
    stream << "radius: " << particle.m_radius << " , mass = " << particle.m_mass << " , color = " <<  particle.m_color
           << " , cell_id = " << particle.m_cell_id << std::endl;
    return stream;
}
#endif
