#ifndef CONVERSIONS_HPP
#define CONVERSIONS_HPP
#include "configuration.hpp"

inline integration_method string_to_integration_method (std::string my_string)
{
    integration_method output;
    if (my_string.compare("euler") == 0)
        output = euler;
    else if (my_string.compare("leap_frog") == 0)
        output = leap_frog;
    else
    {
        std::cerr << "integration method " << my_string << " unkown!" << std::endl;
        throw std::logic_error("unkown integration method");
    }
    return output;
}

inline int string_to_int (std::string my_string)
{
    std::stringstream s_stream(my_string);
    int output;
    s_stream >> output;
    return output;
}

inline double string_to_double (std::string my_string)
{
    std::stringstream s_stream(my_string);
    double output;
    s_stream >> output;
    return output;
}

#endif
