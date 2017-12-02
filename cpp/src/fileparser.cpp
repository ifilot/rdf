//
// RDF (c) by Ivo Filot <i.a.w.filot@tue.nl>
//
// RDF is licensed under a
// Creative Commons Attribution 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by/4.0/>.
//

#include "fileparser.h"

/**
 * @brief      default constructor
 *
 * @param[in]  _filename  path to filename to read geometry from
 */
FileParser::FileParser(const std::string& _filename) :
filename(_filename) {
    this->read_file();
}

/**
 * @brief      Reads a file.
 */
void FileParser::read_file() {
    this->infile.open(this->filename.c_str());

    std::string line; // string to store line in

    size_t dataset_id = 0;

    while(std::getline(infile, line)) {

        // construct holder object
        size_t nr_atoms = boost::lexical_cast<size_t>(line);
        this->datasets.emplace_back(new std::vector<Position>());
        this->datasets.back()->resize(nr_atoms);

        // read contents
        std::getline(infile, line); // discard line
        this->read_atoms(nr_atoms, dataset_id);
        dataset_id++;
    }
}

/**
 * @brief      Reads atoms from file.
 *
 * @param[in]  nr_lines    Number of lines
 * @param[in]  dataset_id  Identifier of the dataset
 */
void FileParser::read_atoms(size_t nr_lines, size_t dataset_id) {

    std::cout << "Reading " << nr_lines << " atoms for dataset #"
              << (dataset_id + 1) << std::endl;

    // construct float_parser object to efficiently read the floats
    // from the lines
    float_parser p;

    std::string line; // string to store line in
    for(size_t i = 0; i < nr_lines; i++) {

        // read the fline
        std::getline(infile, line);

        // set iterators
        std::string::const_iterator b = line.begin();
        std::string::const_iterator e = line.end();

        // parse
        std::vector<float> floats;
        boost::spirit::qi::phrase_parse(b, e, p, boost::spirit::ascii::space, floats);

        // store in the dataset
        this->datasets.back()->at(i).x = floats[1];
        this->datasets.back()->at(i).y = floats[2];
        this->datasets.back()->at(i).z = floats[3];
    }
}
