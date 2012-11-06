#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include    <vector>
#include    <string>

using namespace std;


class Histogram{
    private:
        double  m_start;
        double  m_stop;
        double  m_xToBin;
        unsigned int m_total;
        unsigned int  m_bins;
        vector<unsigned int> m_data;
    public:
        Histogram(const double start, const double stop, const unsigned int bins);

        inline void Clear() {m_data.clear(); m_total = 0;};
        void    Add(const double x);
        void    Dump(string filename);

        inline
        double  GetX(const unsigned int id) {return (static_cast<double>(id)+0.5)/static_cast<double>(m_bins)*(m_stop-m_start);};
        inline
        double  GetXUpper(const unsigned int id) {return (static_cast<double>(id)+1.0)/static_cast<double>(m_bins)*(m_stop-m_start);};
        inline
        unsigned int  GetValue(const unsigned int id) {return m_data[id];};
        inline
        double  GetBinWidth() {return (m_stop-m_start)/static_cast<double>(m_bins);};
};

#endif // HISTOGRAM_H
