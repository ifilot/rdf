// RDF (c) by Ivo Filot <i.a.w.filot@tue.nl>

// RDF is licensed under a
// Creative Commons Attribution 4.0 International License.

// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by/4.0/>.

#ifndef _FLOAT_PARSER_H
#define _FLOAT_PARSER_H

#include <vector>
#include <boost/spirit/include/qi.hpp>

struct float_parser : boost::spirit::qi::grammar<std::string::const_iterator,
                                                  std::vector<float>(),
                                                  boost::spirit::ascii::space_type> {

    float_parser() : float_parser::base_type( vector ) {
        vector  %= +(boost::spirit::qi::float_);
    }

    boost::spirit::qi::rule<std::string::const_iterator, std::vector<float>(), boost::spirit::ascii::space_type> vector;
};

#endif //_FLOAT_PARSER_H
