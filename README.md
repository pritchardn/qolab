# Qolab

Qolab is a C + Intel MKL simulation package/library for the Quantum Approximate Optimisation Algorithm (QAOA)

This package implements the major backbone required for high-performant simulations of the QAOA.
Emphasis is placed on the ability to implement a vast number of high-level parameters (optimisation scheme, 
decomposition etc.)
simply.

## Installation/Setup

After implementing the problem-specific cost function code should be compiled using the provided makefile.

Make sure to adjust the makefile to include the path to your MKL installation root

## Documentation
We use doxygen for our full documentation system. In the `/docs` directory we provide a pdf of the latest documentation 
version. However if you would like to build it yourself please run
```
$ doxygen Doxyfile
```
to build an html and latex version of the documentation.

You will also need the `graphviz` to build the included diagrams

## Usage

Users are only required to implement two main problem dependent features
* The cost function (found in `problem_code.c`)
* The cost function data structure (found in `problem_code.h`)

In this way, the user is able to test the performance of the QAOA without having to formulate an exact circuit 
implementation

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please refer to the currently open issues and project to see what is already coming.

## License
[GPLv3] 

## Acknowledgements

This package currently relies on a number of other libraries. We are greatful for their work.
* Intel MKL https://software.intel.com/en-us/mkl
* Nlopt https://nlopt.readthedocs.io/en/latest/