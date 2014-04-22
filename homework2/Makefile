
# $UWHPSC/codes/fortran/multifile1/Makefile5

OBJECTS = taylor_main.o exptaylor.o
.PHONY: clean

exp_output.txt: taylor_main.exe
	./taylor_main.exe > exp_output.txt

taylor_main.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o taylor_main.exe

%.o : %.f90
	gfortran -c $< 

clean:
	rm -f $(OBJECTS) taylor_main.exe

