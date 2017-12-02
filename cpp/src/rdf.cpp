//
// RDF (c) by Ivo Filot <i.a.w.filot@tue.nl>
//
// RDF is licensed under a
// Creative Commons Attribution 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by/4.0/>.
//

#include "rdf.h"

/**
 * @brief      Default constructor
 *
 * @param[in]  _cutoff   cut-off distance
 * @param[in]  _binsize  bin size
 */
RDF::RDF(float _cutoff, float _binsize) :
cutoff(_cutoff),
binsize(_binsize)
{
    this->nr_bins = (size_t)(cutoff / binsize);
}

/**
 * @brief      construct the RDF from a dataset
 *
 * @param[in]  dataset   pointer to vector of positions (from FileParser)
 * @param[in]  filename  path to output file
 */
void RDF::construct_rdf(const std::vector<Position>* dataset, const std::string& filename) {

    // calculate execution times
    auto start = std::chrono::system_clock::now();
    std::cout << "Construct RDF for " << dataset->size() << " atoms." << std::endl;

    std::vector<float> distances;
    size_t length = dataset->size();

    // calculate RDF (using OpenMP)
    #pragma omp parallel
    {
        std::vector<float> vec_private;
        #pragma omp for schedule(dynamic)
        for(size_t i=0; i<length; i++) {
            for(size_t j=i+1; j<length; j++) {
                const float dx = dataset->at(i).x - dataset->at(j).x;
                const float dy = dataset->at(i).y - dataset->at(j).y;
                const float dz = dataset->at(i).z - dataset->at(j).z;

                const float d = std::sqrt(dx * dx + dy * dy + dz * dz);

                vec_private.push_back(d);
            }
        }

        #pragma omp critical
        distances.insert(distances.end(), vec_private.begin(), vec_private.end());
    }

    // calculate execution time
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Completed RDF calculation in " << elapsed_seconds.count() << " seconds." << std::endl;

    std::cout << "Collecting in " << this->nr_bins << " bins." << std::endl;
    std::vector<size_t> rdf_count(this->nr_bins, 0);
    std::vector<float> rdf(this->nr_bins, 0.0);

    // put all distances into bins
    #pragma omp parallel
    {
        std::vector<size_t> rdf_private(this->nr_bins, 0);
        #pragma omp for nowait
        for(size_t i=0; i<distances.size(); i++) {

            if(distances[i] > this->cutoff) {
                continue;
            }

            size_t container = (size_t)(distances[i] / this->binsize);

            rdf_private[container]++;
        }

        #pragma omp critical
        for(size_t i=0; i<this->nr_bins; i++) {
            rdf_count[i] += rdf_private[i];
        }
    }

    // divide number of interactions by shell volume
    for(size_t i=0; i<this->nr_bins; i++) {
        const float r1 = float(i-1) * this->binsize;                        // inner radius
        const float r2 = float(i) * this->binsize;                          // outer radius
        const float V = 4./3. * M_PI * (std::pow(r2,3) - std::pow(r1,3));   // shell volume
        rdf[i] = (float)rdf_count[i] / V;
    }

    // write results to file
    std::cout << "Writing results to " << filename << std::endl;
    std::ofstream outfile;
    outfile.open(filename);
    for(size_t i=0; i<this->nr_bins; i++) {
        outfile << ((float)i + 0.5) * this->binsize << "\t" << rdf[i] << std::endl;
    }
    outfile.close();
}
