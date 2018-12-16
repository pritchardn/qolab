# Qolab

Qolab is a C + Intel MKL simulation package for the Quantum Approximate Optimisation Algorithm (QAOA)

## Installation


## Usage

Users are only required to implement two main problem dependent features
* The cost function (found in `problem_code.c`)
* The cost function data structure (found in `problem_code.h`)

In this way, the user is able to test the performance of the QAOA without having to formulate an exact circuit implementation

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please refer to the currently open issues and project to see what is already coming.

## License
[GPLv3] 