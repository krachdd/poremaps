OBJECTS= read_problem_params.o geometry.o output.o parallelization.o initialize.o calc_flow.o solver.o boundary_conditions.o evaluation.o
CLINKER=path/to/my/mpiCC

FD_stokes_solver: $(OBJECTS) main_stokes_solver.cc
	-$(CLINKER) -o ../bin/POREMAPS $(OBJECTS) main_stokes_solver.cc -lm
	rm -f $(OBJECTS)

read_problem_params.o: read_problem_params.h read_problem_params.cc
	-$(CLINKER) -o read_problem_params.o -c read_problem_params.cc

geometry.o: constants.h geometry.h geometry.cc
	-$(CLINKER) -o geometry.o -c geometry.cc -lm

output.o: output.h output.cc
	-$(CLINKER) -o output.o -c output.cc -lm

parallelization.o: parallelization.h parallelization.cc
	-$(CLINKER) -o parallelization.o -c parallelization.cc 

initialize.o: initialize.h initialize.cc 
	-$(CLINKER) -o initialize.o -c initialize.cc -lm

calc_flow.o: constants.h solver.h boundary_conditions.h output.h calc_flow.h calc_flow.cc 
	-$(CLINKER) -o calc_flow.o -c calc_flow.cc -lm

solver.o: constants.h solver.h solver.cc
	-$(CLINKER) -o solver.o -c solver.cc 

boundary_conditions.o: constants.h boundary_conditions.h boundary_conditions.cc
	-$(CLINKER) -o boundary_conditions.o -c boundary_conditions.cc

evaluation.o: constants.h evaluation.h evaluation.cc
	-$(CLINKER) -o evaluation.o -c evaluation.cc
