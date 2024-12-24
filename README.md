# IMCCompiler
SIMD IMC compiler that performs design space exploration to find Pareto-optimal designs in terms of netlist size and memory footprint.

## Dependencies
In order to use the compiler, you will need a Windows machine with the following three tools installed:

Mockturtle(https://github.com/lsils/mockturtle)

Gurobi(https://www.gurobi.com/)

Z3(https://github.com/Z3Prover/z3)

Please change the file path in utils.h according to your setting.

We recommend the users to use Visual Studio for this project.

Please add the include directories and library directories of the tools into your project.

Then build the project with Visual Studio.

Note that you may need to use C++17 Language Standard and also update C++ Preprocessor setting in your project to use Mockturle.

Please make sure you can run Mockturle, Gurobi, and Z3 before using our compiler.

## Run
To apply the compiler on an XMG netlist, please use the .v / .bliff / .aig file of the benchmark and put the name of the file in benchmarks.txt (you can optionally add the bound for the netlist).

An example of router benchmark can be found in the files.
