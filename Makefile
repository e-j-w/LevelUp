CFLAGS   = -I./levelup_functions -lm -o2 -Wall

all: levelup

levelup: levelup.c levelup.h
	@echo Making levelup...
	gcc levelup.c $(CFLAGS) -o levelup 
	@echo Tidying up...
	rm -rf *~

clean:
	@echo Cleaning up...
	rm -rf *~ *.o levelup *tmpdatafile*
