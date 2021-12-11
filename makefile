s:	wave_serial.c
	gcc -o wave_serial wave_serial.c -lm

p:	wave_parallel.c
	gcc -std=c99 -o wave_parallel wave_parallel.c -lm -fopenmp

clean:
	rm -f wave_parallel wave_serial rm -f *.o

rs:
	./wave_serial

rp:
	./wave_parallel 8