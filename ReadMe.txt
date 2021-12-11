Instuctions to run code
1) Unzip file contents
2) write make clean to remove .exe files

For serial code
1) write make s
2) run the code using "./wave_serial M T"
	M = number of rows and columns in mesh (single value)
	T = number of time steps
3) write make clean to remove executable

For parallel code
1) write make p
2) run the code using "./wave_parallel N M T"
	N = number of threads to use
	M = number of rows and columns in mesh (single value)
	T = number of time steps
3) write make clean to remove executable

Output data can be processed by using either of the following files:
	wave_serial_plot.m
	wave_parallel_plot.m
