# C++ RDF

If speed is an issue, you can also use the C++ variant of the code.

## Compilation
To compile `rdf`, use the following well-known recipe:

* Make a seperate build directory `mkdir build && cd build`
* Initialize CMake `cmake ../src`
* Build the program `make -j5`

## Usage
To execute, do
```
./rdf -i <path to geometry file> -o <path to output file>
```

## Visualization
You can visualize the result, for instance, using Gnuplot. In Gnuplot, use the following:
```
plot "output.file" u 1:2 w l
```
