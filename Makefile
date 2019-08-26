all:
	g++ main.cpp -o main

eqs:
	g++ Equation.cpp -c -o Equation.o
	g++ Mathieu.cpp -c -o Mathieu.o