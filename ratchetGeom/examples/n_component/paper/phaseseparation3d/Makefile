numcomp = 6
ROOT_DIR =	../../../..
INCS	=	-I$(ROOT_DIR)/src
DEFS    =	-g -O3 -Wno-deprecated -Wall -std=c++17 -fopenmp -DNUMCOMPONENTS=$(numcomp) -DTHREEDIMENSIONAL #-DMPIPARALLEL#-DTHREEDIMENSIONAL 
LIBS	=	-lm
COMP	=	mpicxx
DEPS    =   main.hh
OBJ     =	main.cc
strip_templates = 1
ifdef strip_templates
STRIP_TEMPLATES = -fdiagnostics-color=always 2>&1 | sed 's/ \[with [^]]*\]//g'
endif

all: run$(numcomp).exe 

%.o: %.cc $(DEPS)
	$(COMP) -c -o $@ $< $(DEFS) $(STRIP_TEMPLATES) $(LIBS)

run$(numcomp).exe: $(OBJ)
	$(COMP) -o run$(numcomp).exe $(OBJ) $(INCS) $(DEFS) $(STRIP_TEMPLATES) $(LIBS)

clean:
	find . -type f -name \*.exe -exec rm {} \;
	

#rm -r data