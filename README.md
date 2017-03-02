# Genetic
basic genetic algorithm implemented in C to calculate maximum point for a given function on a given function in [1,10] at the accuracy of 4 decimal places
-- updated on 2017.03.02
--to get the maximum point, call the gene_algorithm(int mode). if mode == 0, it will return the maximum point; it will return the minimum point if mode != 0. then you can call func(double value), and it will return the function value of the given value.
if you want to calculate another interval other than [1,10], you need to update upper_bound and lower_bound and let chromosome_size satisfy this formula:
2^(chromosome_size-1) < (upper_bound - lower_bound) * 10 ^ accuracy < 2^chromosome_size
if you want to modify accuracy, just update chromosome_size;
if you want to modify the function, replace the function body of func();
if you want to modify the size of population, update popu_size;

good luck.
