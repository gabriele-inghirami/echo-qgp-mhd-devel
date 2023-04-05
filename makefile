# *** ECHO-QGP makefile ***

# compiler to be called for producing a "serial" executable (i.e. to be run on a single cpu)
com_ser=gfortran

# compiler to be called for producing a "parallel" executable (i.e. to be run on multi core cpus) using the MPI library
com_par=mpif90

# flag specific to the gfortran compiler to use 8 bytes precision as the default for "real" types
flags=-fdefault-real-8 -fdefault-double-8

# optimization or debugging flags

# safe optimization flag, valid for all processors
opt=-O3

# optimization flag that enables specific cpu features that can speed up the computations
# CAREFUL: this option is safe only if the processors of the computer where the code has been compiled and of the computer where it will run are equal or equivalent 
# so, this option is safe when the compilation and the execution of the code happen on the same computer, but it may be not safe, for example, in a cluster
opt+=-march=native

# optimization flag for amd barcelona processors 
#opt=-O3 -march=barcelona

# debugging flags for GNU gdb (and also ddd, Affinic Debugger, Allinea DDT,...)
#opt=-O0 -g -ggdb #-fbounds-check

# debugging flag for Absoft FX3 debugger
#opt=-g -gdwarf-3

#additional debugging flags to detect floating point exceptions
#opt+=-ffpe-trap=invalid,zero,overflow

#additional flag for gfortran >= 10
flags+=-fallow-argument-mismatch

# standard default compilation target: the serial version of the executable
all:	  echo.exe 

obj_ser = parallel_nompi.o holib_s.o eos_s.o common_s.o viscous_s.o system_s.o work_s.o evolve_s.o hypersurface_s.o out_s.o glaubermc_s.o init_s.o echo_s.o

echo.exe:	  ${obj_ser}
		  ${com_ser} -o echo.exe ${obj_ser} ${flags} ${opt}

parallel_nompi.o:	  parallel_nompi.f90 makefile
		  ${com_ser} ${flags} ${opt} -c parallel_nompi.f90 -o parallel_nompi.o
holib_s.o:	  holib.f90 makefile
		  ${com_ser} ${flags} ${opt} -c holib.f90 -o holib_s.o
eos_s.o:	  eos.f90 makefile
		  ${com_ser} ${flags} ${opt} -c eos.f90 -o eos_s.o
common_s.o:	  common.f90 makefile
		  ${com_ser} ${flags} ${opt} -c common.f90 -o common_s.o
viscous_s.o:	  viscous.f90 makefile
		  ${com_ser} ${flags} ${opt} -c viscous.f90 -o viscous_s.o
system_s.o:	  system.f90 makefile
		  ${com_ser} ${flags} ${opt} -c system.f90 -o system_s.o
work_s.o:	  work.f90 common.f90 holib.f90 system.f90 makefile
		  ${com_ser} ${flags} ${opt} -c work.f90 -o work_s.o
evolve_s.o:	  evolve.f90 parallel_nompi.f90 common.f90 work.f90 system.f90 makefile
		  ${com_ser} ${flags} ${opt} -c evolve.f90 -o evolve_s.o
hypersurface_s.o:	  hypersurface.f90 parallel_nompi.f90 common.f90 makefile
		  ${com_ser} ${flags} ${opt} -c hypersurface.f90 -o hypersurface_s.o
out_s.o:	  out.f90 parallel_nompi.f90 common.f90 eos.f90 makefile
		  ${com_ser} ${flags} ${opt} -c out.f90 -o out_s.o
glaubermc_s.o:	  glaubermc.f90 parallel_nompi.f90 common.f90 system.f90 makefile
		  ${com_ser} ${flags} ${opt} -c glaubermc.f90 -o glaubermc_s.o
init_s.o:	  init.f90 parallel_nompi.f90 common.f90 system.f90 glaubermc.f90 makefile
		  ${com_ser} ${flags} ${opt} -c init.f90 -o init_s.o
echo_s.o:	  echo.f90 parallel_nompi.f90 common.f90 evolve.f90 out.f90 makefile
		  ${com_ser} ${flags} ${opt} -c echo.f90 -o echo_s.o

# target: the parallel (MPI) version of the executable
par:	  	  echo_par.exe 

obj_par = parallel_mpi.o holib_p.o eos_p.o common_p.o viscous_p.o system_p.o work_p.o evolve_p.o hypersurface_p.o out_p.o glaubermc_p.o init_p.o echo_p.o

echo_par.exe:	  ${obj_par}
		  ${com_par} -o echo.exe ${flags} ${obj_par} ${opt}


parallel_mpi.o:	  parallel_mpi.f90 makefile
		  ${com_par} ${flags} ${opt} -c parallel_mpi.f90 -o parallel_mpi.o
holib_p.o:	  holib.f90 makefile
		  ${com_par} ${flags} ${opt} -c holib.f90 -o holib_p.o
eos_p.o:	  eos.f90 makefile
		  ${com_par} ${flags} ${opt} -c eos.f90 -o eos_p.o
common_p.o:	  common.f90 makefile
		  ${com_par} ${flags} ${opt} -c common.f90 -o common_p.o
viscous_p.o:	  viscous.f90 makefile
		  ${com_par} ${flags} ${opt} -c viscous.f90 -o viscous_p.o
system_p.o:	  system.f90 makefile
		  ${com_par} ${flags} ${opt} -c system.f90 -o system_p.o
work_p.o:	  work.f90 common.f90 holib.f90 system.f90 makefile
		  ${com_par} ${flags} ${opt} -c work.f90 -o work_p.o
evolve_p.o:	  evolve.f90 parallel_mpi.f90 common.f90 work.f90 system.f90 makefile
		  ${com_par} ${flags} ${opt} -c evolve.f90 -o evolve_p.o
hypersurface_p.o:	  hypersurface.f90 parallel_mpi.f90 common.f90 makefile
		  ${com_par} ${flags} ${opt} -c hypersurface.f90 -o hypersurface_p.o
out_p.o:	  out.f90 parallel_mpi.f90 common.f90 eos.f90 makefile
		  ${com_par} ${flags} ${opt} -c out.f90 -o out_p.o
glaubermc_p.o:	  parallel_mpi.f90 common.f90 system.f90 makefile
		  ${com_par} ${flags} ${opt} -c glaubermc.f90 -o glaubermc_p.o
init_p.o:	  init.f90 parallel_mpi.f90 common.f90 system.f90 glaubermc.f90 makefile
		  ${com_par} ${flags} ${opt} -c init.f90 -o init_p.o
echo_p.o:	  echo.f90 parallel_mpi.f90 common.f90 evolve.f90 out.f90 makefile
		  ${com_par} ${flags} ${opt} -c echo.f90 -o echo_p.o

# target: some postprocessing tools to analyze the results 
tools:	  tools.exe

obj_tools = parallel_nompi.o holib_s.o eos_s.o

tools.exe:	  ${obj_tools}
		  ${com_ser} ${flags} ${opt} tools/timev.f90 -o tools/timev.exe
		  ${com_ser} ${flags} ${opt} tools/fromecho.f90 -o tools/fromecho.exe
		  ${com_ser} ${flags} ${opt} tools/fromecho2d.f90 -o tools/fromecho2d.exe
		  ${com_ser} ${flags} ${opt} tools/readx.f90 -o tools/readx.exe


# target: to delete the products of the compilation (executable, object files, modules...)
clean:
		  \rm *.mod *.o *.exe

# target: to delete both the files produced by the compiler and the files produced by an echo-qgp run
cleanall:
		  \rm -rf *.mod *.o *.exe log wound.dat thick.dat thick_params.dat ed.dat summary.dat *.err *.out out/* graph/* postproc/readx/* outr* partcoll.dat tools/*.exe
