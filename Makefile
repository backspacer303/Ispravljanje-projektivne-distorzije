PROGRAM = projekcija
CC      = g++
FLAGS   = -lm -lpthread -lX11 -ljpeg -lpng

$(PROGRAM): projekcija.cpp
	$(CC) projekcija.cpp -o $(PROGRAM) $(FLAGS)
