//
// RDF (c) by Ivo Filot <i.a.w.filot@tue.nl>
//
// RDF is licensed under a
// Creative Commons Attribution 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by/4.0/>.
//

#ifndef _FILE_PARSER_H
#define _FILE_PARSER_H

#include <string>
#include <fstream>
#include <iostream>

#include <boost/lexical_cast.hpp>

#include "float_parser.h"
#include "position.h"

/**
 * @brief      Parses a geometry file and stores it in an internal object
 */
class FileParser {

private:
    std::string filename;
    std::ifstream infile;
    std::vector<std::unique_ptr<std::vector<Position>>> datasets;

public:

    /**
     * @brief      default constructor
     *
     * @param[in]  _filename  path to filename to read geometry from
     */
    FileParser(const std::string& _filename);

    inline const std::vector<Position>* get_dataset(int i) const {
        if(i < 0) {
            return this->datasets.back().get();
        }

        if(i >= this->datasets.size()) {
            exit(-1);
        }

        return this->datasets[i].get();
    }

private:

    /**
     * @brief      Reads a file.
     */
    void read_file();

    /**
     * @brief      Reads atoms from file.
     *
     * @param[in]  nr_lines    Number of lines
     * @param[in]  dataset_id  Identifier of the dataset
     */
    void read_atoms(size_t nr_lines, size_t dataset_id);

};

#endif //_FILE_PARSER_H
