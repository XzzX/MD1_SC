#include "plot_utils.hpp"


void drawCellSubdivision(const unsigned int n_horizontal, const unsigned int n_vertical, const double width, const double height, const cColor &color, const Window &w)
{
    cColor red = {1.0, 0.0, 0.0};

    for (unsigned int i = 0; i < n_horizontal; i++)
        w.drawHline((static_cast<double>(i)/static_cast<double>(n_horizontal) - 1.0 / 2.0) * width, height, red);

    for (unsigned int i = 0; i < n_vertical; i++)
        w.drawVline(width, (static_cast<double>(i)/static_cast<double>(n_vertical) - 1.0 / 2.0) * height, red);
}

void drawValues(const std::pair<std::list<double>::const_iterator, std::list<double>::const_iterator> &values, arrangement &pos, const cColor& color, const Window &w)
{
    if (values.first == values.second)
        return;

    double y_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::min();
    double x_min = 0;
    double x_max = 0;;
    if (pos.fixed == false)
    {
        //search for minimum and maximum values

        unsigned int counter = 0;
        for (std::list<double>::const_iterator it = values.first; it != values.second; it++)
        {
            if (*it < y_min)
                y_min = *it;
            else if (*it > y_max)
                y_max = *it;
            counter++;
        }
        x_max = counter;

        pos.x_range_min = 0.0;
        pos.x_range_max = x_max;

        pos.y_range_min = y_min;
        pos.y_range_max = y_max;
        pos.fixed = true;
    } else {
        y_min = pos.y_range_min;
        y_max = pos.y_range_max;
        x_min = pos.x_range_min;
        x_max = pos.x_range_max;
    }

    double delta_x = 50.0/(x_max - x_min);
    double delta_x_over_4 = delta_x / 4.0;

    unsigned index = 0;
    for (std::list<double>::const_iterator it = values.first; it != values.second; it++)
    {
        Vector v1(delta_x * index - 25.0, (*it - y_min)*35.0/(y_max-y_min) - 20.0, 0.0);
        w.drawCross(v1, delta_x_over_4, color);
        index++;
    }
}

void drawCoordinateAxis(arrangement &pos, const cColor& color, const Window &w)
{
    //w.drawArrow(Vector(pos.x_range_min), Vector(), color);
    //w.drawArrow(Vector(), Vector(), color);
}
