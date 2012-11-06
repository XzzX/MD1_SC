#ifndef WINDOW_H
#define WINDOW_H
#include <SFML/Window.hpp>
#include <stdexcept>
#include <GL/glut.h>

#include <iostream>
#include <list>
#include "Vector.hpp"
#include "particles.hpp"

class Window
{
  public:
    Window();
    Window(std::string title, Vector _L, const unsigned int subdivisions);
    ~Window();
    bool isOpen();
    void processEvents();

    void drawCircle(const Vector &pos, const cColor &color, const double r) const;
    void drawRectangle(const Vector &pos, const double width, const double height, const cColor& color) const;
    void drawHline(const double x, const double y, const cColor &color) const;
    void drawVline(const double x, const double y, const cColor &color) const;
    void drawCross(const Vector &x, const double size, const cColor &color) const;
    void drawErrorbar(const Vector &x, const double size, const double error, const cColor &color) const;
    void drawVector(const Vector &vec, const cColor &color, const double zoom) const;
    void drawArrow(const Vector &x, const Vector &y, const cColor &color) const;
    void display();
    void close();
    void clear();

    int q_value;
    int w_value;
    int a_value;
    int s_value;
    int y_value;
    int x_value;

	void SetCameraPosition(const Vector& pos);

  private:
    GLuint circleList;
    Vector L;
    Vector position_camera;
    Vector position_mouse;

    double windowRatio;
    sf::Window* App;
    unsigned int referenceCounter;
    unsigned int m_segments;
};

inline Window::Window()
{
  App=NULL;
}


inline double dmod(double a, double b)
{
  return a-static_cast<int>(a/b)*b;
}

inline void Window::display()
{
    if (!App)
        throw std::logic_error("Window::isOpen(): there is no window!");
        // Finally, display rendered frame on screen
    App->Display();
}

inline void Window::close()
{
    App->Close();
}

inline Window::~Window()
{
    if (App)
    {
        // FIXME: copy-constructor not dealing right with this:
        delete App;
    }
}

inline bool Window::isOpen()
{
    if (!App)
        throw std::logic_error("Window::isOpen(): there is no window!");
    return App->IsOpened();
}


#endif
