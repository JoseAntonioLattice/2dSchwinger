FC = gfortran
TARGET = 2dU1.exe

SRC = src
BIN = bin

SOURCE = indices.f90 statistics.f90 pbc.f90 configurations.f90 arrays.f90 parameters.f90 starts.f90 u1.f90 Gradient_flow.f90 dirac.f90 CG.f90 hmc.f90 observables.f90 dynamics.f90 main.f90

OBJECT = $(patsubst %, $(BIN)/%, $(SOURCE:.f90=.o ) )

#FFLAGS = -Wall -Wextra -fcheck=all -O0 -J$(BIN) -I$(BIN) -cpp
FFLAGS = -ffree-line-length-512 -O3 -J$(BIN) -I$(BIN) -cpp #-DPARALLEL 



$(BIN)/$(TARGET): $(OBJECT)
	$(FC) -o $@ $^ -L ~/Fortran/lib -lnum2str -lfiles -lconstants

$(BIN)/%.o: $(SRC)/%.f90
	$(FC) -I ~/Fortran/include $(FFLAGS) -c $< -o $@ -L ~/Fortran/lib -lnum2str  -lfiles

.PHONY: help run clean

run:
#	@echo "input/input_parameters.nml" | time $(BIN)/$(TARGET)
#	{ echo "input/input_parameters.nml"; echo 2 2; } | cafrun -n 4 $(BIN)/$(TARGET)
	@echo "input/input_parameters.nml" | LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/Fortran/lib time $(BIN)/$(TARGET)
#	{ echo "input/input_parameters.nml"; echo 4 1; } | LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/Fortran/lib cafrun -n 4 $(BIN)/$(TARGET)

clean:
	rm -f $(OBJECT) $(BIN)/$(TARGET)

help:
	@echo "src: $(SOURCE)"
