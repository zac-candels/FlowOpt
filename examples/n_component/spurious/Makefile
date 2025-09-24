ROOT_DIR =	../../..
INCS	=	-I$(ROOT_DIR)/src
DEFS    =	-g -O3 -Wno-deprecated -Wall -std=c++17 -fopenmp #-DMPIPARALLEL #-DTHREEDIMENSIONAL #
LIBS	=	-lm
COMP	=	mpicxx
DEPS    =       
OBJ     =	main.cc
strip_templates = 1
ifdef strip_templates
STRIP_TEMPLATES = -fdiagnostics-color=always 2>&1 | sed 's/ \[with [^]]*\]//g'
endif

all: run.exe 

%.o: %.cc $(DEPS)
	$(COMP) -c -o $@ $< $(DEFS) $(STRIP_TEMPLATES) $(LIBS)

run.exe: $(OBJ)
	$(COMP) -o run.exe $(OBJ) $(INCS) $(DEFS) $(STRIP_TEMPLATES) $(LIBS)

clean:
	rm -f *.o run.exe
	

#rm -r data