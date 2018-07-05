# Ising lattices complete enumeration of dilution effects

This set of scripts performs a set of calculations to get the dependence of residual entropy on field for kagome, triangular and pyrochlore lattices.

There are three folders with propitiate files: `kagome`, `triangular` and `pyrochlore`.
And three `*.cpp` scripts inside. 
Scripts asks two parameters: number of holes (zero spins, which are magnetically neutral) and field. These are command prompt parameters, not command-line options! Just run the program in console and follow prompts on screen.

Then, the program enumerate all possible variants of dilution (permutations), enumerate all possible configurations of a set of Ising spins for each permutation, build DOS (density of states), and write the number of ground states for each permutation to the file.

The resulting file is named `g_%1_%2.dat`, where `%1` - previously entered number of holes, `%2` - external field.

## Running scripts
There are two scripts placed along with cpp file: `run.sh` and `res_average.sh`.

For using it, firstly compile your program to file named `main.o`

Example: `g++ kagome.cpp -O3 -o main.o`

`run.sh` runs main.o program with varied parameters of field and dilution. Inside the file, `holes=( 0 3 7 10 15 )` variable is responsible for variants of possible number of holes in sample. Variable `h=$(LANG=en_US seq 0.0 0.01 5.0)` is the enumeration of possible field values. In the example it goes from 0.0 to 5.0 with 0.1 step. `LANG=en_US` is responsible for correct comma/dot decimal separator (bash-specific feature).
The result of calculation will be in `res` folder.

`res_average.sh` makes the averagement of calculated data, and saves the result to file `g_aver_%1.dat` where `%1` is the number of holes. After, you can plot the data. There are two same variables `holes` and `h` inside of `res_average.sh` file like in `run.sh`. You can tune them.