#include "particles.hpp"

cColor createColor(unsigned int i, unsigned int N)
{
    cColor output;
    const double H=i*360/N;
    const int hi=static_cast<int>(H/60.0);
    const double f=H/60.0-hi;
    const double S=1.0;
    const double V=1; // 255 for full range
    const double p=V*(1-S);
    const double q=V*(1-S*f);
    const double t=V*(1-S*(1-f));
    if (hi==0 || hi==6)
        compoment_to_color(V,t,p,output);
    else if (hi==1)
        compoment_to_color(q,V,p,output);
    else if (hi==2)
        compoment_to_color(p,V,t,output);
    else if (hi==3)
        compoment_to_color(p,q,V,output);
    else if (hi==4)
        compoment_to_color(t,p,V,output);
    else if (hi==5)
        compoment_to_color(V,p,q,output);
    else
        compoment_to_color(1.0,1.0,1.0,output);
    return output;
}

