CFLAGS   = -I./levelup_functions -lm -o2 -Wall -lreadline

all: levelup

levelup: levelup.c levelup.h levelup_functions/*
	@echo Making levelup...
	gcc levelup.c $(CFLAGS) -o levelup 
	@echo Tidying up...
	rm -rf *~

clean:
	@echo Cleaning up...
	rm -rf *~ *.o levelup *tmpdatafile*
