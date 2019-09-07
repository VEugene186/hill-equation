all: tools
	g++ --std=c++11 --openmp Equation.o Mathieu.o Mathieu2.o LiouvilleHill.o EulerPoisson.o RungeKutta.o EulerPoissonFixedPoints.o EulerPoissonTracker.o main.cpp -o main

eqs:
	g++ --std=c++11 Equation.cpp -c -o Equation.o
	g++ --std=c++11 Mathieu.cpp -c -o Mathieu.o
	g++ --std=c++11 Mathieu2.cpp -c -o Mathieu2.o
	g++ --std=c++11 LiouvilleHill.cpp -c -o LiouvilleHill.o
	g++ --std=c++11 EulerPoisson.cpp -c -o EulerPoisson.o

method: eqs
	g++ --std=c++11 RungeKutta.cpp -c -o RungeKutta.o

tools: method
	g++ --std=c++11 EulerPoissonFixedPoints.cpp -c -o EulerPoissonFixedPoints.o
	g++ --std=c++11 EulerPoissonTracker.cpp -c -o EulerPoissonTracker.o
clean:
	rm -rf *.o main