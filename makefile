big_bang_simulator:	main.o universe.o
	g++ -o big_bang_simulator main.o universe.o -lGL -lGLU -lglut

main.o: main.cpp universe.h
	g++ -c -Wall main.cpp

universe.cpp: universe.h
	g++ -c -Wall universe.cpp
