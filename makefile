all: compilar

compilar: silly_season.cpp
	g++ -O2 -o silly_season silly_season.cpp

clean: 
	rm silly_season