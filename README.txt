This GitHub-repository contains all the source code for the numerical library
IGA_Soliton. It can be used for solving several partial differential equations
used for describing the motion of solitary waves. The code is written entirely
in C++, using some specialized libraries. The graph plotting is implemented with
Python. Further details on using the software can be found in the instructions
manual which is included in this repository.


--------------------------------------------------------------------------------
| Package requirements                                                         |
--------------------------------------------------------------------------------

To compile and run the source code, the following packages are required:
  * C++ programming language (for computation)
    - GNU Compiler Collection 6.3.0
    - MPI 1.10
    - OpenMP 4.5
    - PETSc 3.7
    - TinyXML 2.6.2
    - Expreval
  * Python programming language (for plotting graphs)
    - NumPy
    - Matplotlib
  * Cmake 3.5.0 (cross-platform compilation)

PETSc and TinyXML are libraries which must be installed on the computer, while
Expreval only needs to be included in the repository itself.

The easiest way of including NumPy and Matplotlib is to install the Anaconda
launcher, a complete program containing all the relevant Python-libraries for
scientific computation.

https://sourceforge.net/projects/tinyxml/
https://sourceforge.net/projects/expreval/

--------------------------------------------------------------------------------
| Compiling and running procedure                                              |
--------------------------------------------------------------------------------

To create the executable files in each folder, use the shellscript-commands:
    mkdir build && cd build
    CC=gcc CXX=g++ cmake ..
    make
