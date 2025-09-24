ROOT_DIR =	../../..
INCS	=	-I$(ROOT_DIR)/src
DEFS    =	-g -O3 -Wno-deprecated -Wall -std=c++17 -fopenmp #-DMPIPARALLEL#
LIBS	=	-lm
COMP	=	mpicxx
DEPS    =       
OBJ     =	main.cc

define strip_templates
endef
ifdef strip_templates
STRIP_TEMPLATES = -fdiagnostics-color=always 2>&1 | sed 's/ \[with [^]]*\]//g'
endif

all: run.exe 

two: run2.exe

kinetic: runkinetic.exe

%.o: %.cc $(DEPS)
	$(COMP) -c -o $@ $< $(DEFS) $(LIBS)

run.exe: $(OBJ)
	$(COMP) -o run.exe $(OBJ) $(INCS) $(DEFS) $(LIBS) $(STRIP_TEMPLATES)

run2.exe: $(OBJ)
	$(COMP) -o run2.exe $(OBJ) $(INCS) $(DEFS) $(LIBS) $(STRIP_TEMPLATES)

runkinetic.exe: $(OBJ)
	$(COMP) -o runkinetic.exe $(OBJ) $(INCS) $(DEFS) -DKINETIC $(LIBS) $(STRIP_TEMPLATES)

clean:
	rm -f *.o run.exe

cleankinetic:
	rm -f *.o runkinetic.exe

#rm -r data