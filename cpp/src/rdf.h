//
// RDF (c) by Ivo Filot <i.a.w.filot@tue.nl>
//
// RDF is licensed under a
// Creative Commons Attribution 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by/4.0/>.
//

#ifndef _RDF_H
#define _RDF_H

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

#include "position.h"

/**
 * @brief      Class for calculating Radial Distribution Function
 */
class RDF {

private:
    float cutoff;       // cutoff distance
    float binsize;      // bin size
    size_t nr_bins;     // number of bins

public:

    /**
     * @brief      Default constructor
     *
     * @param[in]  _cutoff   cut-off distance
     * @param[in]  _binsize  bin size
     */
    RDF(float _cutoff, float _binsize);

    /**
     * @brief      construct the RDF from a dataset
     *
     * @param[in]  dataset   pointer to vector of positions (from FileParser)
     * @param[in]  filename  path to output file
     */
    void construct_rdf(const std::vector<Position>* dataset, const std::string& filename);

private:

};

#endif //_RDF_H
