# RDF

Calculates the Radial Distribution Function (RDF) given a unit cell and atoms inside the unit cell.

Currently hardcoded unit cells:
* FCC
* BCC
* HCP
* SC

If speed is an issue, use the [C++ variant](cpp). The C++ variant uses OpenMP. On a i7-4790K, the Python version needs about 4 minutes to parse 13,824 atoms. The C++ variant only uses 1.5 seconds ()
