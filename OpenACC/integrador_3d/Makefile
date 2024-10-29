FC = gfortran
LIBS = -llapack -lblas
BIN = integrador
OBJECTS = main.o sp_integrador_2d.o get_coef.o
FFLAGS1 = -o
FFLAGS2 = -c

$(BIN) : $(OBJECTS)
	$(FC) $(FFLAGS1) $(BIN) $(OBJECTS)

mod_constantes.o :
	$(FC) $(FFLAGS2) mod_constantes.f90

main.o : main.f90 mod_constantes.o
	$(FC) $(FFLAGS2) main.f90

sp_integrador_2d.o : sp_integrador_2d.f90 mod_constantes.o
	$(FC) $(FFLAGS2) sp_integrador_2d.f90

get_coef.o : get_coef.f90 mod_constantes.o
	$(FC) $(FFLAGS2) get_coef.f90

.PHONY : clean
clean : 
	rm $(BIN) *.o *.mod
