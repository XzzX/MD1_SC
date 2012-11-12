#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include    <vector>
#include    <string>

using namespace std;


class Histogram{
    private:
        ///left border of the histogram
        double  m_start;
        ///right border of the histogram
        double  m_stop;
        ///number of bins
        unsigned int    m_bins;
        ///conversion factor between x value and corresponding bin number
        double  m_xToBin;
        ///total amount of entries in histogram
        unsigned int    m_total;
        ///amount of events per bin
        vector<unsigned int> m_data;
    public:
        Histogram(const double start, const double stop, const unsigned int bins);

        ///Clear all the data of the histogram
        inline void Clear() {m_data.clear(); m_total = 0;};
        ///Add a new Event
        void    Add(const double x);
        ///Save the histogram to a file
        void    Dump(string filename);

        ///Get middle position of the bin
        inline
        double  GetX(const unsigned int id) {return (static_cast<double>(id)+0.5)/static_cast<double>(m_bins)*(m_stop-m_start);};
        ///get right position of the bin
        inline
        double  GetXUpper(const unsigned int id) {return (static_cast<double>(id)+1.0)/static_cast<double>(m_bins)*(m_stop-m_start);};
        ///get amount of events in bin
        inline
        unsigned int  GetValue(const unsigned int id) {return m_data[id];};
        ///get width of a bin
        inline
        double  GetBinWidth() {return (m_stop-m_start)/static_cast<double>(m_bins);};
};

#endif // HISTOGRAM_H
