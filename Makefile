all: method
	g++ --openmp Equation.o Mathieu.o RungeKutta.o main.cpp -o main

eqs:
	g++ Equation.cpp -c -o Equation.o
	g++ Mathieu.cpp -c -o Mathieu.o

method: eqs
	g++ RungeKutta.cpp -c -o RungeKutta.o

clean:
	rm -rf *.o main