#include "Histogram.h"

#include    <fstream>

Histogram::Histogram(const double start, const double stop, const unsigned int bins):
    m_start(start), m_stop(stop), m_bins(bins), m_xToBin(static_cast<double>(bins)/(stop-start)), m_total(0){
    m_data.resize(bins);
}

void    Histogram::Add(const double x){
    unsigned int i = (x-m_start) * m_xToBin;
    if (i>=m_bins) return ;
    m_data[i]++;
    m_total++;
}

void    Histogram::Dump(string filename){
    fstream fout(filename.c_str(), fstream::out);

    for(unsigned int i=0; i<m_bins; i++)
        fout << (static_cast<double>(i)+0.5)/static_cast<double>(m_bins)*(m_stop-m_start) << "\t" << m_data[i] << std::endl;

    fout.close();
}
