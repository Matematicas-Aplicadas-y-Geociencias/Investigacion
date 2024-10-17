FC = gfortran
LIBS = -llapack -lblas
BIN = integrador
OBJECTS = main.o sp_integrador_2d.o resolver_sel.o
FFLAGS1 = -o
FFLAGS2 = -c

$(BIN) : $(OBJECTS)
	$(FC) $(FFLAGS1) $(BIN) $(OBJECTS) $(LIBS)

mod_constantes.o :
	$(FC) $(FFLAGS2) mod_constantes.f90

main.o : main.f90 mod_constantes.o
	$(FC) $(FFLAGS2) main.f90

sp_integrador_2d.o : sp_integrador_2d.f90 mod_constantes.o
	$(FC) $(FFLAGS2) sp_integrador_2d.f90

resolver_sel.o : resolver_sel.f90 mod_constantes.o
	$(FC) $(FFLAGS2) resolver_sel.f90

.PHONY : clean
clean : 
	rm $(BIN) *.o *.mod
