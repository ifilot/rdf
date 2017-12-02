//
// RDF (c) by Ivo Filot <i.a.w.filot@tue.nl>
//
// RDF is licensed under a
// Creative Commons Attribution 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by/4.0/>.
//

#include <iostream>
#include <tclap/CmdLine.h>

#include "fileparser.h"
#include "rdf.h"

int main(int argc, char* argv[]) {

    try {

        TCLAP::CmdLine cmd("Calculate RDF from file.", ' ', "1.0");

        // input file
        TCLAP::ValueArg<std::string> arg_input_filename("i","input","Input file (i.e. geom.dat)",true,"__NONE__","filename");
        cmd.add(arg_input_filename);

        // output file
        TCLAP::ValueArg<std::string> arg_output_filename("o","output","Output file (i.e. rdf.dat)",true,"__NONE__","filename");
        cmd.add(arg_output_filename);

        cmd.parse(argc, argv);

        std::string input_filename = arg_input_filename.getValue();
        std::string output_filename = arg_output_filename.getValue();

        FileParser fp(input_filename);

        RDF rdf(15, 0.05);
        rdf.construct_rdf(fp.get_dataset(-1), output_filename);

    } catch (TCLAP::ArgException &e) {

        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;

    }

    return 0;
}
