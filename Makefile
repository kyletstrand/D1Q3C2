CFLAGS = -g -std=c99 -I $(HOME)/include 
LDFLAGS = -g -L $(HOME)/lib -L /usr/X11R6/lib 
LIBS = -lgraph -lm -lX11

% : %.c

% : %.o 
	$(CC) $(LDFLAGS) $< $(LIBS) -o $@

%.o : %.c
	$(CC) $(CFLAGS) -c $<  -o $@