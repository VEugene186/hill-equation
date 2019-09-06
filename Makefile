all: method
	g++ --std=c++11 --openmp Equation.o Mathieu.o Mathieu2.o LiouvilleHill2.o EulerPoisson.o RungeKutta.o main.cpp -o main

eqs:
	g++ --std=c++11 Equation.cpp -c -o Equation.o
	g++ --std=c++11 Mathieu.cpp -c -o Mathieu.o
	g++ --std=c++11 Mathieu2.cpp -c -o Mathieu2.o
	g++ --std=c++11 LiouvilleHill2.cpp -c -o LiouvilleHill2.o
	g++ --std=c++11 EulerPoisson.cpp -c -o EulerPoisson.o

method: eqs
	g++ --std=c++11 RungeKutta.cpp -c -o RungeKutta.o

clean:
	rm -rf *.o main