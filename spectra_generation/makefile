all: thermal 
# merge particles


# serial or openmp (choose and uncomment only one of the two)
# serial mode
#com = gfortran -O3 -fdefault-real-8 -march=native
# with openmp
com = gfortran -O3 -fdefault-real-8 -march=native -fopenmp

# serial for debugging purposes
#com = gfortran -O0 -g -ggdb3 -fno-strict-overflow -fdefault-real-8 -fcheck=bounds

obj_th = bessel.o common.o init.o io.o thermal.o
		
thermal:	${obj_th}
		${com} ${obj_th} -o thermal.exe	

bessel.o:		bessel.f90 makefile
			${com} -c bessel.f90
common.o:		common.f90  makefile
			${com} -c common.f90
init.o:	 		init.f90 common.f90 makefile
			${com} -c init.f90 
			
io.o:			io.f90 init.f90 common.f90  makefile
			${com} -c io.f90			

thermal.o:		thermal.f90 common.f90 init.f90 io.f90 makefile
			${com} -c thermal.f90
			

clean:		
		\rm *.mod *.o  thermal.exe 

cleanall:	
		\rm *.mod *.o *D.exe log *~ *.err *.out thermal.exe


		
