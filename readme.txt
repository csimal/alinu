# readme file for conjugate gradient in C

# Author
csimal

# Required Libraries ########################################
Due to CLAPACK being deprecated in favor of LAPACKE, we use the
latter in this project. To install LAPACKE on Linux Mint/Ubuntu,
simply run:

sudo apt install liblapacke-dev

# Compilation ###############################################
To compile the project, simply run "make"
Use "make clean" to delete every compiled file.

# Execution #################################################
The makefile generates an executable named "main"
This program does all the things required for the analysis of the
various methods, including:
* solving all three systems using the conjugate gradient method
first without, then with various preconditionners, and outputting
the number of iterations needed to converge on the screen, and 
writing the iterates to files in the same directory

* computing the eigenvalues of every matrix used and writing them to
a file in the same directory
