! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 * 
! *                                                                           *         
! *  Version: 1.5.0-apha                                                      *
! *                                                                           *
! *  Copyright (C) 2015-2019 The ECHO-QGP team                                * 
! *                                                                           *
! *  File: out.f90                                                            *
! *                                                                           *
! *  License: GPL version 2.0 (Please, read the file LICENSE.TXT)             *
! *                                                                           *
! *  This program is free software; you can redistribute it and/or            *
! *  modify it under the terms of the GNU General Public License              *
! *  as published by the Free Software Foundation; either version 2           *
! *  of the License, or (at your option) any later version.                   *
! *                                                                           *
! *  This program is distributed in the hope that it will be useful,          *
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
! *  GNU General Public License for more details.                             *
! *                                                                           *
! *  You should have received a copy of the GNU General Public License        *
! *  along with this program; if not, write to the Free Software              *
! *  Foundation, Inc., 51 Franklin Street, Fifth Floor,                       *
! *  Boston, MA  02110-1301, USA.                                             *
! *                                                                           *
! *  Authors: Luca Del Zanna (luca.delzanna@unifi.it)                         *
! *           Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)            *
! *                                                                           *
! *  Contributors:                                                            *
! *                                                                           *
! *  Acknowledgments:                                                         *
! *                                                                           *
! *****************************************************************************

! development note: careful to the mpi_realtype if we switch from/to 4/8 bytes real precision

module out


  use parallel
  use common
  use eos
  use hypersurface
    
  implicit none

  integer ldir
  integer itime1,itime2,iratio
  integer output
  real(KIND=8), allocatable, dimension(:,:,:) :: arr_to_print_8
  real(KIND=8), allocatable, dimension(:,:,:) :: arr_to_print_8_bis
  real(KIND=8), allocatable, dimension(:,:,:) :: arr_to_print_8_ter
  real(KIND=8), allocatable, dimension(:,:,:) :: arr_to_print_8_quater
  real(KIND=8), allocatable, dimension(:,:,:) :: arr_to_print_8_quinquies

contains

! *****************************************************************************

subroutine out_alloc
!-- Allocate in/out arrays
  implicit none
  integer :: allocate_result
  integer :: i !just for counting
  
  if (pe0) then
    if (id_of_the_run<10000) write (dir_out,'(a4,i4,a1)') 'outr' ,id_of_the_run,'/'
    if (id_of_the_run<1000) write (dir_out,'(a5,i3,a1)') 'outr0' ,id_of_the_run,'/'
    if (id_of_the_run<100) write (dir_out,'(a6,i2,a1)') 'outr00' ,id_of_the_run,'/'
    if (id_of_the_run<10) write (dir_out,'(a7,i1,a1)') 'outr000', id_of_the_run,'/'
    commandline_outdir_length=len_trim(adjustl(custom_output_directory))
    if(commandline_outdir_length .gt. 0) then
      allocate(character(len=len_trim(adjustl(custom_output_directory))+10) :: prefix_dir)
      prefix_dir=trim(adjustl(custom_output_directory))//"/"//dir_out
    else
      allocate(character(len=9) :: prefix_dir)
      prefix_dir=dir_out
    end if

    if(out_freeze) then
      call allocate_hypersurface_arrays
    end if

    allocate(arr_to_print_8(nx,ny,nz), arr_to_print_8_bis(nx,ny,nz), arr_to_print_8_ter(nx,ny,nz), stat=allocate_result)
    if(allocate_result /=0) then
      write(*,*)  "Proc.0 - Error, I can't allocate the arr_to_print_8 arrays into out_alloc"
      write(*,*)  "(Problem in source file out.f90, subroutine out_alloc)"
      call exit(1)
    end if
    if(viscosity) then
      allocate(arr_to_print_8_quater(nx,ny,nz), arr_to_print_8_quinquies(nx,ny,nz), stat=allocate_result)
      if(allocate_result /=0) then
        write(*,*)  "Proc.0 - Error, I can't allocate the arr_to_print_8_quater/quinquies array into out_alloc"
        write(*,*)  "(Problem in source file out.f90, subroutine out_alloc)"
        call exit(1)
      end if
    end if
  end if

  !even if these arrays are used only in parallel runs, since they are small for simplicity we always allocate them 
  allocate(recvcounts(0:npe-1), displs(0:npe-1), stat=allocate_result)
  if(allocate_result /=0) then
    write(*,*) "Error, I can't allocate recvcounts or displs array for the MPI_Gatherv funcion (source file common.f08)"
    call exit(1)
  end if

  if((restart_type .ne. 0) .and. (pe0)) then
    allocate(uall(nx,ny,nz,1:nv), vallold(nx,ny,nz,krh:kvold_end), uallold(nx,ny,nz,krh:kvold_end), stat=allocate_result)
    if(allocate_result /=0) then
      write(*,*)  "Proc.0 - Error, I can't allocate 'uall', 'vallold' or/and 'uallold' arrays into out_alloc"
      write(*,*)  "(Problem in source file out.f90, subroutine out_alloc)"
      call exit(1)
    end if
    uall=0.
    vallold=0.
    uallold=0. 
  end if

 displs(0)=0
 if (prl) then
  do i=0,npe-1
    if (i .lt. ipec) then
      recvcounts(i)=nx/npe+1
    else
      recvcounts(i)=nx/npe
    end if
    if(i .gt. 0) displs(i)=displs(i-1)+recvcounts(i-1)
  end do
 end if
 
 if(pe0 .and. derivatives_out) then
   allocate(derivatives_all(nx,ny,nz,1:nderivatives), stat=allocate_result)
      if(allocate_result /=0) then
      write(*,*)  "Proc.0 - Error, I can't allocate 'derivatives_all' array into out_alloc"
      write(*,*)  "(source file out.f08)"
      call exit(1)
      end if
      
 end if 

end subroutine out_alloc

! *****************************************************************************

subroutine out_vars_summary()
!-- Write data on file 'vars_summary.dat'
  use viscous
       
  if (pe0) then
    if (id_of_the_run<10000) write (dir_out,'(a4,i4,a1)') 'outr' ,id_of_the_run,'/'
    if (id_of_the_run<1000) write (dir_out,'(a5,i3,a1)') 'outr0' ,id_of_the_run,'/'
    if (id_of_the_run<100) write (dir_out,'(a6,i2,a1)') 'outr00' ,id_of_the_run,'/'
    if (id_of_the_run<10) write (dir_out,'(a7,i1,a1)') 'outr000', id_of_the_run,'/'
    if(commandline_outdir_length .gt. 0) then
      prefix_dir=trim(adjustl(custom_output_directory))//"/"//dir_out
    else
      prefix_dir=dir_out
    end if
  
    open  (13,file=prefix_dir//'vars_summary.dat',form='formatted',status='replace')
        013     format(i1,a35)
        write(13,*) 'Summary of variable printed in the output files. 1 means printed, 0 not printed'
        if(output_precision .eq. 4) then
          write(13,013) 4, 'precision'
        else
          write(13,013) 8, 'precision'
        end if
        write(13,013) out_sel%density, 'density'
        write(13,013) out_sel%vx,'vx'
        write(13,013) out_sel%vy,'vy'
        write(13,013) out_sel%vz,'vz'
        write(13,013) out_sel%pressure,'pressure'
        write(13,013) out_sel%energy_density, 'energy_density'
        write(13,013) out_sel%temperature, 'temperature'
        write(13,013) out_sel%entropy_density, 'entropy_density'
        if (viscosity) then
          write(13,013) out_sel%bulk, 'bulk viscosity'
          write(13,013) out_sel%pitt, 'pi^tt shear tensor component'
          write(13,013) out_sel%pitx, 'pi^tx shear tensor component'
          write(13,013) out_sel%pity, 'pi^ty shear tensor component'
          write(13,013) out_sel%pitz, 'pi^tz shear tensor component'
          write(13,013) out_sel%pixy, 'pi^xy shear tensor component'
          write(13,013) out_sel%pixz, 'pi^xz shear tensor component'
          write(13,013) out_sel%piyz, 'pi^yz shear tensor component'
          write(13,013) out_sel%pixx, 'pi^xx shear tensor component'
          write(13,013) out_sel%piyy, 'pi^yy shear tensor component'
          write(13,013) out_sel%pizz, 'pi^zz shear tensor component'
        else
          write(13,013) 0, 'bulk viscosity'
          write(13,013) 0, 'pi^tt shear tensor component'
          write(13,013) 0, 'pi^tx shear tensor component'
          write(13,013) 0, 'pi^ty shear tensor component'
          write(13,013) 0, 'pi^tz shear tensor component'
          write(13,013) 0, 'pi^xy shear tensor component'
          write(13,013) 0, 'pi^xz shear tensor component'
          write(13,013) 0, 'pi^yz shear tensor component'
          write(13,013) 0, 'pi^xx shear tensor component'
          write(13,013) 0, 'pi^yy shear tensor component'
          write(13,013) 0, 'pi^zz shear tensor component'
        end if
        write(13,013) out_sel%v0, 'vt (or gamma Lorentz factor)'
        if (mhd) then
          write(13,013) out_sel%bx, 'Bx'
          write(13,013) out_sel%by, 'By'
          write(13,013) out_sel%bz, 'Bz'
          write(13,013) out_sel%ex, 'Ex'
          write(13,013) out_sel%ey, 'Ey'
          write(13,013) out_sel%ez, 'Ez'
          if(divclean) write(13,013) out_sel%glm, 'GLM'
          write(13,013) out_sel%rc, 'el. charge density in comoving frame'
        else
          write(13,013) 0, 'Bx'
          write(13,013) 0, 'By'
          write(13,013) 0, 'Bz'
          write(13,013) 0, 'Ex'
          write(13,013) 0, 'Ey'
          write(13,013) 0, 'Ez'
          write(13,013) 0, 'GLM'
          write(13,013) 0, 'el. charge density in comoving frame'
        end if

     close (13)
  
  end if

end subroutine out_vars_summary

! *****************************************************************************

subroutine out_vars(iout,nooverwrite)
!-- Write data on file 'out000n.dat'
  use system_eqgp, only: g_cov,system_EXB,system_compute_E_in_ideal_MHD
  use eos
  use evolve
  use viscous
  integer,intent(inout) :: iout
  integer      :: lv,i,lvm,lvp,proc
  integer      :: iii,jjj,kkk,lll,mmm,nnn
  real(8) energy_density_value
  real(8) temperature_value
  real(8) entropy_density_value
  real(8),dimension(1:3) :: E_field
  character(len=7) :: fname
  integer :: nooverwrite

  logical :: retrieve_derived_data
  integer :: error_flag

077 format(e14.8)

  error_flag=0
  if((out_sel%energy_density .eq. 1) .or. (out_sel%temperature .eq. 1) .or. (out_sel%entropy_density .eq. 1)) then
    retrieve_derived_data=.true.
  else
    retrieve_derived_data=.false.
  end if

  if (pe0) then
    if (iout<10000) write (fname,'(a3,i4)') 'out' ,iout
    if (iout<1000) write (fname,'(a4,i3)') 'out0' ,iout
    if (iout<100) write (fname,'(a5,i2)') 'out00' ,iout
    if (iout<10) write (fname,'(a6,i1)') 'out000', iout
    if (id_of_the_run<10000) write (dir_out,'(a4,i4,a1)') 'outr' ,id_of_the_run,'/'
    if (id_of_the_run<1000) write (dir_out,'(a5,i3,a1)') 'outr0' ,id_of_the_run,'/'
    if (id_of_the_run<100) write (dir_out,'(a6,i2,a1)') 'outr00' ,id_of_the_run,'/'
    if (id_of_the_run<10) write (dir_out,'(a7,i1,a1)') 'outr000', id_of_the_run,'/'
    if(commandline_outdir_length .gt. 0) then
      prefix_dir=trim(adjustl(custom_output_directory))//"/"//dir_out
    else
      prefix_dir=dir_out
    end if
       
    open  (10,file=prefix_dir//fname//'.dat',form='unformatted', access='stream',status='replace')

    if(nooverwrite .eq. 0) then
       open (21,file=prefix_dir//'time.dat',form='formatted',status='replace')
    else
       open (21,file=prefix_dir//'time.dat',form='formatted',status='old',access='append')
    end if

    call system_clock(itime1)
    if(output_precision .eq. 8) then
       write (10) t
    else
       write (10) real(t,4)
    end if
    
    write (21,"(I5,f16.9)") iout,t
    close (21)

  end if

   
  if((out_sel%rc .eq. 1) .and. mhd) call compute_rho_comoving_frame
  
  call collect_distributed_arrays

  if(pe0) then
   if(output_precision .eq. 8) then
     if(out_sel%density .eq. 1) then
       write(10) vall(:,:,:,krh)
     end if
     if(out_sel%vx .eq. 1) then
       write(10) vall(:,:,:,kvx)
     end if
     if(out_sel%vy .eq. 1) then
       write(10) vall(:,:,:,kvy)
     end if
     if(out_sel%vz .eq. 1) then
       write(10) vall(:,:,:,kvz)
     end if
     if(out_sel%pressure .eq. 1) then
       write(10) vall(:,:,:,kpr)
     end if
     if(retrieve_derived_data) then
       do nnn=1,nz
        do mmm=1,ny
         do lll=1,nx
          call get_derived_data(vall(lll,mmm,nnn,krh),vall(lll,mmm,nnn,kpr), energy_density_value,&
          &temperature_value, entropy_density_value,error_flag )
          if(error_flag .gt. 0) then
            write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
            errcode=error_flag
            write(*,*) "Error code:", errcode
            write(*,*) "Position on the grid: ix=",lll,"iy=",mmm,"iz=",nnn
            run_crashed=.true.
            return
          end if
          arr_to_print_8(lll,mmm,nnn)=energy_density_value
          arr_to_print_8_bis(lll,mmm,nnn)=temperature_value
          arr_to_print_8_ter(lll,mmm,nnn)=entropy_density_value
         end do
        end do
       end do
       if(out_sel%energy_density .eq. 1) write(10) arr_to_print_8
       if(out_sel%temperature .eq. 1) write(10) arr_to_print_8_bis
       if(out_sel%entropy_density .eq. 1) write(10) arr_to_print_8_ter
      end if !end if retrieve data
      if(viscosity) then
       if (out_sel%bulk .eq. 1) then
         write(10) vall(:,:,:,kpibu)
       end if
       if((out_sel%pitx .eq. 1) .or. (out_sel%pity .eq. 1) .or. (out_sel%pitt .eq. 1) .or. &
         & (out_sel%pitz .eq. 1) ) then
         if(obtained .eq. 'no') then
           do nnn=1,nz
            do mmm=1,ny
             do lll=1,nx
              call get_derived_pi(vall(lll,mmm,nnn,:),arr_to_print_8(lll,mmm,nnn),arr_to_print_8_bis(lll,mmm,nnn),&
              &arr_to_print_8_ter(lll,mmm,nnn),arr_to_print_8_quater(lll,mmm,nnn))
             end do
            end do
           end do
         else
           do nnn=1,nz
            do mmm=1,ny
             do lll=1,nx
              call get_derived_pi_zz(vall(lll,mmm,nnn,:),arr_to_print_8(lll,mmm,nnn),arr_to_print_8_bis(lll,mmm,nnn),&
              &arr_to_print_8_ter(lll,mmm,nnn),arr_to_print_8_quater(lll,mmm,nnn),arr_to_print_8_quinquies(lll,mmm,nnn))
             end do
            end do
           end do
         end if
        end if
        if(out_sel%pitt .eq. 1) write(10) arr_to_print_8
        if(out_sel%pitx .eq. 1) write(10) arr_to_print_8_bis
        if(out_sel%pity .eq. 1) write(10) arr_to_print_8_ter
        if(out_sel%pitz .eq. 1) write(10) arr_to_print_8_quater
        if(out_sel%pixy .eq. 1) then
           write(10) vall(:,:,:,kpixy)
        end if
        if(out_sel%pixz .eq. 1) then
           write(10) vall(:,:,:,kpixz)
        end if
        if(out_sel%piyz .eq. 1) then
           write(10) vall(:,:,:,kpiyz)
        end if
        if(out_sel%pixx .eq. 1) then
           write(10) vall(:,:,:,kpixx)
        end if
        if(out_sel%piyy .eq. 1) then
           write(10) vall(:,:,:,kpiyy)
        end if
        if(out_sel%pizz .eq. 1) then
           write(10) arr_to_print_8_quinquies
        end if
      end if !end if viscosity
      if(out_sel%v0 .eq. 1) then
        arr_to_print_8=1./sqrt(1.-vall(:,:,:,kvx)**2.*g_cov(1)-vall(:,:,:,kvy)**2.*g_cov(2)-vall(:,:,:,kvz)**2.*g_cov(3))
        write(10) arr_to_print_8
      end if
      if(mhd) then
        if(out_sel%bx .eq. 1) then
           write(10) vall(:,:,:,kbx)
        end if
        if(out_sel%by .eq. 1) then
           write(10) vall(:,:,:,kby)
        end if
        if(out_sel%bz .eq. 1) then
           write(10) vall(:,:,:,kbz)
        end if
        if(rmhd) then
           if(out_sel%ex .eq. 1) then
             write(10) vall(:,:,:,kex)
           end if
           if(out_sel%ey .eq. 1) then
             write(10) vall(:,:,:,key)
           end if
           if(out_sel%ez .eq. 1) then
             write(10) vall(:,:,:,kez)
           end if
        else
         if((out_sel%ex .eq. 1) .or. (out_sel%ey .eq. 1) .or. (out_sel%ez .eq. 1)) then
           do nnn=1,nz
             do mmm=1,ny
             do lll=1,nx
              call system_compute_E_in_ideal_MHD(vall(lll,mmm,nnn,kvx:kvz), vall(lll,mmm,nnn,kbx:kbz), E_field)
              arr_to_print_8(lll,mmm,nnn)=E_field(1)
              arr_to_print_8_bis(lll,mmm,nnn)=E_field(2)
              arr_to_print_8_ter(lll,mmm,nnn)=E_field(3)/g_cov(3) !because the previous function returns E covariant
             end do
            end do
           end do
           if(out_sel%ex .eq. 1) write(10) arr_to_print_8
           if(out_sel%ey .eq. 1) write(10) arr_to_print_8_bis
           if(out_sel%ez .eq. 1) write(10) arr_to_print_8_ter
         end if !end check any E field component to print
        end if!end else rmhd 
        if(divclean .and. (out_sel%glm .eq. 1)) then
          write(10) vall(:,:,:,kglm)
        end if
        if(out_sel%rc .eq. 1) then
          write(10) vall(:,:,:,krc)
        end if
      end if !end if mhd
     close (10)
    else!output precision 4 bytes
     if(out_sel%density .eq. 1) then
       write(10) real(vall(:,:,:,krh),4)
     end if
     if(out_sel%vx .eq. 1) then
       write(10) real(vall(:,:,:,kvx),4)
     end if
     if(out_sel%vy .eq. 1) then
       write(10) real(vall(:,:,:,kvy),4)
     end if
     if(out_sel%vz .eq. 1) then
       write(10) real(vall(:,:,:,kvz),4)
     end if
     if(out_sel%pressure .eq. 1) then
       write(10) real(vall(:,:,:,kpr),4)
     end if
     if(retrieve_derived_data) then
       do nnn=1,nz
        do mmm=1,ny
         do lll=1,nx
          call get_derived_data(vall(lll,mmm,nnn,krh),vall(lll,mmm,nnn,kpr), energy_density_value,&
          &temperature_value, entropy_density_value,error_flag )
          if(error_flag .gt. 0) then
            write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
            errcode=error_flag
            write(*,*) "Error code:", errcode
            write(*,*) "Position on the grid: ix=",lll,"iy=",mmm,"iz=",nnn
            run_crashed=.true.
            return
          end if
          arr_to_print_8(lll,mmm,nnn)=energy_density_value
          arr_to_print_8_bis(lll,mmm,nnn)=temperature_value
          arr_to_print_8_ter(lll,mmm,nnn)=entropy_density_value
         end do
        end do
       end do
       if(out_sel%energy_density .eq. 1) write(10) real(arr_to_print_8,4)
       if(out_sel%temperature .eq. 1) write(10) real(arr_to_print_8_bis,4)
       if(out_sel%entropy_density .eq. 1) write(10) real(arr_to_print_8_ter,4)
      end if !end if retrieve data
      if(viscosity) then
       if (out_sel%bulk .eq. 1) then
         write(10) real(vall(:,:,:,kpibu),4)
       end if
       if((out_sel%pitx .eq. 1) .or. (out_sel%pity .eq. 1) .or. (out_sel%pitt .eq. 1) .or. &
         & (out_sel%pitz .eq. 1) ) then
         if(obtained .eq. 'no') then
           do nnn=1,nz
            do mmm=1,ny
             do lll=1,nx
              call get_derived_pi(vall(lll,mmm,nnn,:),arr_to_print_8(lll,mmm,nnn),arr_to_print_8_bis(lll,mmm,nnn),&
              &arr_to_print_8_ter(lll,mmm,nnn),arr_to_print_8_quater(lll,mmm,nnn))
             end do
            end do
           end do
         else
           do nnn=1,nz
            do mmm=1,ny
             do lll=1,nx
              call get_derived_pi_zz(vall(lll,mmm,nnn,:),arr_to_print_8(lll,mmm,nnn),arr_to_print_8_bis(lll,mmm,nnn),&
              &arr_to_print_8_ter(lll,mmm,nnn),arr_to_print_8_quater(lll,mmm,nnn),arr_to_print_8_quinquies(lll,mmm,nnn))
             end do
            end do
           end do
         end if
        end if
        if(out_sel%pitt .eq. 1) write(10) real(arr_to_print_8,4)
        if(out_sel%pitx .eq. 1) write(10) real(arr_to_print_8_bis,4)
        if(out_sel%pity .eq. 1) write(10) real(arr_to_print_8_ter,4)
        if(out_sel%pitz .eq. 1) write(10) real(arr_to_print_8_quater,4)
        if(out_sel%pixy .eq. 1) then
           write(10) real(vall(:,:,:,kpixy),4)
        end if
        if(out_sel%pixz .eq. 1) then
           write(10) real(vall(:,:,:,kpixz),4)
        end if
        if(out_sel%piyz .eq. 1) then
           write(10) real(vall(:,:,:,kpiyz),4)
        end if
        if(out_sel%pixx .eq. 1) then
           write(10) real(vall(:,:,:,kpixx),4)
        end if
        if(out_sel%piyy .eq. 1) then
           write(10) real(vall(:,:,:,kpiyy),4)
        end if
        if(out_sel%pizz .eq. 1) then
           write(10) real(arr_to_print_8_quinquies,4)
        end if
      end if !end if viscosity
      if(out_sel%v0 .eq. 1) then
        arr_to_print_8=1./sqrt(1.-vall(:,:,:,kvx)**2.*g_cov(1)-vall(:,:,:,kvy)**2.*g_cov(2)-vall(:,:,:,kvz)**2.*g_cov(3))
        write(10) real(arr_to_print_8,4)
      end if
      if(mhd) then
        if(out_sel%bx .eq. 1) then
           write(10) real(vall(:,:,:,kbx),4)
        end if
        if(out_sel%by .eq. 1) then
           write(10) real(vall(:,:,:,kby),4)
        end if
        if(out_sel%bz .eq. 1) then
           write(10) real(vall(:,:,:,kbz),4)
        end if
        if(rmhd) then
          if(out_sel%ex .eq. 1) then
             write(10) real(vall(:,:,:,kex),4)
          end if
          if(out_sel%ey .eq. 1) then
             write(10) real(vall(:,:,:,key),4)
          end if
          if(out_sel%ez .eq. 1) then
             write(10) real(vall(:,:,:,kez),4)
          end if
        else
          if((out_sel%ex .eq. 1) .or. (out_sel%ey .eq. 1) .or. (out_sel%ez .eq. 1)) then
           do nnn=1,nz
            do mmm=1,ny
             do lll=1,nx
              call system_compute_E_in_ideal_MHD(vall(lll,mmm,nnn,kvx:kvz), vall(lll,mmm,nnn,kbx:kbz), E_field)
              arr_to_print_8(lll,mmm,nnn)=E_field(1)
              arr_to_print_8_bis(lll,mmm,nnn)=E_field(2)
              arr_to_print_8_ter(lll,mmm,nnn)=E_field(3)/g_cov(3) !because the previous function returns E covariant
             end do
            end do
           end do
           if(out_sel%ex .eq. 1) write(10) real(arr_to_print_8,4)
           if(out_sel%ey .eq. 1) write(10) real(arr_to_print_8_bis,4)
           if(out_sel%ez .eq. 1) write(10) real(arr_to_print_8_ter,4)
          end if!end out_sel%ex...
        end if !end if rmhd
        if(divclean .and. (out_sel%glm .eq. 1 )) then
          write(10) real(vall(:,:,:,kglm),4)
        end if
        if(out_sel%rc .eq. 1) then
          write(10) real(vall(:,:,:,krc),4)
        end if
      end if !end if mhd
     close (10)
    end if !end if precision 4
  endif !end if pe0

  if(derivatives_out) call out_derivatives(iout,1)

  if(pe0) then
    if(flows_out) call compute_on_transverse_plane(iout,error_flag)
  end if

  if(pe0) then
    call system_clock(itime2,iratio)
    print '(3x,2(a,f12.8),a)','Time:',t,' - '//fname//'.dat',real((itime2-itime1),8)/real(iratio,8),' secs'
  end if 

end subroutine out_vars
! *****************************************************************************

subroutine compute_rho_comoving_frame
use evolve
implicit none


integer i, j, k, il, ih, jl, jh, kl, kh
real(KIND=8) :: time_now, time_past, dt, dx, dy, dz, g
real(KIND=8) :: Ex_now, Ey_now, Ez_now, Ex_past, Ey_past, Ez_past, Ex_now_p, Ex_now_m, Ey_now_p, Ey_now_m, Ez_now_p, Ez_now_m
real(KIND=8) :: dBxcov_dy, dBxcov_dz, dBycov_dx, dBycov_dz, dBzcov_dx, dBzcov_dy, dExcont_dx
real(KIND=8) :: dEycont_dy, dEzcont_dz, dExcont_dt, dEycont_dt, dEzcont_dt, Jx, Jy, Jz, divE, curlBx, curlBy, curlBz, Jdotv
! here we define the covariant four velocity, its derivatives and the right/left cell four vectors to compute the derivatives
real(KIND=8) :: uu0,uu1,uu2,uu3,du01,du02,du03,du12,du13,du21,du23,du31,du32,du10,du20,du30
real(KIND=8) :: u0rx,u0ry,u0rz,u1ry,u1rz,u2rx,u2rz,u3rx,u3ry!values of u covariant components at the right cell along x,y or z
real(KIND=8) :: u0lx,u0ly,u0lz,u1ly,u1lz,u2lx,u2lz,u3lx,u3ly!values of u covariant components at the left cell along x,y or z
real(KIND=8) :: g_past,uu0_past,uu1_past,uu2_past,uu3_past!values of the gamma Lorentz factor and 4 velocity components in the past
real(KIND=8) :: grx,glx,gry,gly,grz,glz!values of the gamma Lorentz factor at right and left cell along x,y or z
logical, save :: not_first_time=.false.
logical, parameter :: maxwell=.false.
real(KIND=8), parameter :: e0123=1.d0, e1023=-1.d0, e0132=-1.d0, e0312=1.d0, e3012=-1.d0, e3021=1.d0, e3201=-1.d0, e3210=1.d0
real(KIND=8), parameter :: e0321=-1.d0, e0231=1.d0, e0213=-1.d0, e2013=1.d0, e2031=-1.d0, e1032=1.d0, e1320=1.d0, e1302=-1.d0
real(KIND=8), parameter :: e1230=-1.d0, e1203=1.d0, e2310=-1.d0, e2301=1.d0, e3120=-1.d0, e3102=1.d0, e2130=1.d0, e2103=-1.d0
real(KIND=8) :: om0,om1,om2,om3!components of the kinematic vorticity four-vector
real(KIND=8) :: vdotB,bb0,bb1,bb2,bb3!v_k B^k and components of the magnetic field four-vector as in eq. 36 of arXiv:1609.03042


!if the time has not changed, we don't recompute rho comov
if(t==time_comp_rho_comov) return
if(t==tmin) return

if(coordinates .eq. 1) then
  time_now=1
  time_past=1
else
  time_now=t
  time_past=timeold
end if

dt=t-timeold

dz=2.d0/ddz(1)
dx=2.d0/ddx(1)
dy=2.d0/ddy(1)

if((nx==1) .or. (ny==1)) then
  write(*,*) "Sorry, but the computation of the charge density works only with both nx and ny > 1"
  call exit(4)
end if 

w(ix1:ix2,iy1:iy2,iz1:iz2,1:kvold_end)=v(ix1:ix2,iy1:iy2,iz1:iz2,1:kvold_end)
call evolve_bcx(ix1,iy1,iz1,1,kvold_end,kv(1:kvold_end))
w(ix1:ix2,iy1:iy2,iz1:iz2,1:kvold_end)=v(ix1:ix2,iy1:iy2,iz1:iz2,1:kvold_end)
call evolve_bcy(ix1,iy1,iz1,1,kvold_end,kv(1:kvold_end))
if (nz>1) then
    w(ix1:ix2,iy1:iy2,iz1:iz2,1:kvold_end)=v(ix1:ix2,iy1:iy2,iz1:iz2,1:kvold_end)
    call evolve_bcz(ix1,iy1,iz1,1,kvold_end,kv(1:kvold_end))

 if(maxwell) then
  do k=iz1,iz2
    do j=iy1,iy2
      do i=ix1,ix2
        dBxcov_dy=(w(i,j+1,k,kbx)-w(i,j-1,k,kbx))/dy
        dBxcov_dz=(w(i,j,k+1,kbx)-w(i,j,k-1,kbx))/dz
        dBycov_dx=(w(i+1,j,k,kby)-w(i-1,j,k,kby))/dx
        dBycov_dz=(w(i,j,k+1,kby)-w(i,j,k-1,kby))/dz
        dBzcov_dx=g_cov(3)*(w(i+1,j,k,kbz)-w(i-1,j,k,kbz))/dx
        dBzcov_dy=g_cov(3)*(w(i,j+1,k,kbz)-w(i,j-1,k,kbz))/dy
        Ex_now=-time_now*(v(i,j,k,kvy)*v(i,j,k,kbz)-v(i,j,k,kvz)*v(i,j,k,kby))
        Ex_now_p=-time_now*(w(i+1,j,k,kvy)*w(i+1,j,k,kbz)-w(i+1,j,k,kvz)*w(i+1,j,k,kby))
        Ex_now_m=-time_now*(w(i-1,j,k,kvy)*w(i-1,j,k,kbz)-w(i-1,j,k,kvz)*w(i-1,j,k,kby))
        Ey_now=-time_now*(v(i,j,k,kvz)*v(i,j,k,kbx)-v(i,j,k,kvx)*v(i,j,k,kbz))
        Ey_now_p=-time_now*(w(i,j+1,k,kvz)*w(i,j+1,k,kbx)-w(i,j+1,k,kvx)*w(i,j+1,k,kbz))
        Ey_now_m=-time_now*(w(i,j-1,k,kvz)*w(i,j-1,k,kbx)-w(i,j-1,k,kvx)*w(i,j-1,k,kbz))
        Ez_now=-(v(i,j,k,kvx)*v(i,j,k,kby)-v(i,j,k,kvy)*v(i,j,k,kbx))/time_now
        Ez_now_p=-(w(i,j,k+1,kvx)*w(i,j,k+1,kby)-w(i,j,k+1,kvy)*w(i,j,k+1,kbx))/time_now
        Ez_now_m=-(w(i,j,k-1,kvx)*w(i,j,k-1,kby)-w(i,j,k-1,kvy)*w(i,j,k-1,kbx))/time_now
        Ex_past=-time_past*(vold(i,j,k,kvy)*vold(i,j,k,kbz)-vold(i,j,k,kvz)*vold(i,j,k,kby))
        Ey_past=-time_past*(vold(i,j,k,kvz)*vold(i,j,k,kbx)-vold(i,j,k,kvx)*vold(i,j,k,kbz))
        Ez_past=-(vold(i,j,k,kvx)*vold(i,j,k,kby)-vold(i,j,k,kvy)*vold(i,j,k,kbx))/time_past
        dExcont_dx=(Ex_now_p-Ex_now_m)/dx
        dEycont_dy=(Ey_now_p-Ey_now_m)/dy
        dEzcont_dz=(Ez_now_p-Ez_now_m)/dz
        dExcont_dt=(time_now*Ex_now-time_past*Ex_past)/dt
        dEycont_dt=(time_now*Ey_now-time_past*Ey_past)/dt
        dEzcont_dt=(time_now*Ez_now-time_past*Ez_past)/dt
        g=1/sqrt(1-v(i,j,k,kvx)**2-v(i,j,k,kvy)**2-v(i,j,k,kvz)**2*g_cov(3))
        divE=dExcont_dx+dEycont_dy+dEzcont_dz
        curlBx=dBzcov_dy-dBycov_dz
        curlBy=dBxcov_dz-dBzcov_dx
        curlBz=dBycov_dx-dBxcov_dy
        Jx=(curlBx-dExcont_dt)/time_now
        Jy=(curlBy-dEycont_dt)/time_now
        Jz=(curlBz-dEzcont_dt)/time_now
        Jdotv=Jx*v(i,j,k,kvx)+Jy*v(i,j,k,kvy)+Jz*v(i,j,k,kvz)*g_cov(3)
        v(i,j,k,krc)=g*(divE-Jdotv)
       end do   
    end do
  end do
 else !we compute the charge density using the vorticity
  do k=iz1,iz2
    do j=iy1,iy2
      do i=ix1,ix2
        g=1.d0/sqrt(1.d0-w(i,j,k,kvx)**2-w(i,j,k,kvy)**2-w(i,j,k,kvz)**2*time_now**2)
        uu0=-g!u0 covariant
        uu1=g*w(i,j,k,kvx)
        uu2=g*w(i,j,k,kvy)
        uu3=g*w(i,j,k,kvz)*time_now**2
        g_past=1.d0/sqrt(1.d0-vold(i,j,k,kvx)**2-vold(i,j,k,kvy)**2-vold(i,j,k,kvz)**2*time_past**2)
        uu0_past=-g_past!u0 covariant at previous timestep
        uu1_past=g_past*vold(i,j,k,kvx)
        uu2_past=g_past*vold(i,j,k,kvy)
        uu3_past=g_past*vold(i,j,k,kvz)*time_past**2
        grx=1.d0/sqrt(1.d0-w(i+1,j,k,kvx)**2-w(i+1,j,k,kvy)**2-w(i+1,j,k,kvz)**2*time_now**2)
        glx=1.d0/sqrt(1.d0-w(i-1,j,k,kvx)**2-w(i-1,j,k,kvy)**2-w(i-1,j,k,kvz)**2*time_now**2)
        gry=1.d0/sqrt(1.d0-w(i,j+1,k,kvx)**2-w(i,j+1,k,kvy)**2-w(i,j+1,k,kvz)**2*time_now**2)
        gly=1.d0/sqrt(1.d0-w(i,j-1,k,kvx)**2-w(i,j-1,k,kvy)**2-w(i,j-1,k,kvz)**2*time_now**2)
        grz=1.d0/sqrt(1.d0-w(i,j,k+1,kvx)**2-w(i,j,k+1,kvy)**2-w(i,j,k+1,kvz)**2*time_now**2)
        glz=1.d0/sqrt(1.d0-w(i,j,k-1,kvx)**2-w(i,j,k-1,kvy)**2-w(i,j,k-1,kvz)**2*time_now**2)
        u1ry=gry*w(i,j+1,k,kvx)
        u1ly=gly*w(i,j-1,k,kvx)
        u1rz=grz*w(i,j,k+1,kvx)
        u1lz=glz*w(i,j,k-1,kvx)
        u2rx=grx*w(i+1,j,k,kvy)
        u2lx=glx*w(i-1,j,k,kvy)
        u2rz=grz*w(i,j,k+1,kvy)
        u2lz=glz*w(i,j,k-1,kvy)
        u3rx=grx*w(i+1,j,k,kvz)*time_now**2
        u3lx=glx*w(i-1,j,k,kvz)*time_now**2
        u3ry=gry*w(i,j+1,k,kvz)*time_now**2
        u3ly=gly*w(i,j-1,k,kvz)*time_now**2
        !we recall that dx, dy and dz are already two times the cell spacing, see their definiton above
        du10=(uu1-uu1_past)/dt!du_x/dt
        du20=(uu2-uu2_past)/dt!du_y/dt
        du30=(uu3-uu3_past)/dt!du_z/dt
        du01=-(grx-glx)/dx!du_0/dx=-d_gammaLorentz/dx
        du02=-(gry-gly)/dy!du_0/dy=-d_gammaLorentz/dy
        du03=-(grz-glz)/dz!du_0/dz=-d_gammaLorentz/dz
        du12=(u1ry-u1ly)/dy!du_x/dy
        du13=(u1rz-u1lz)/dz!du_x/dz
        du21=(u2rx-u2lx)/dx!du_y/dx
        du23=(u2rz-u2lz)/dz!du_y/dz
        du31=(u3rx-u3lx)/dx!du_z/dx
        du32=(u3ry-u3ly)/dy!du_z/dy
        !we use eq. A1 in https://arxiv.org/pdf/1904.01530.pdf
        om0=(e1203*du21*uu3+e2103*du12*uu3+e3102*du13*uu2+e1302*du31*uu2+e2301*du32*uu1+e3201*du23*uu1)/time_now
        om1=(e0213*du20*uu3+e2013*du02*uu3+e0312*du30*uu2+e3012*du03*uu2+e2310*du32*uu0+e3210*du23*uu0)/time_now
        om2=(e0123*du10*uu3+e1023*du01*uu3+e0321*du30*uu1+e3021*du03*uu1+e1320*du31*uu0+e3120*du13*uu0)/time_now
        om3=(e0231*du20*uu1+e2031*du02*uu1+e0132*du10*uu2+e1032*du01*uu2+e1230*du21*uu0+e2130*du12*uu0)/time_now
        vdotB=w(i,j,k,kvx)*w(i,j,k,kbx)+w(i,j,k,kvy)*w(i,j,k,kby)+w(i,j,k,kvz)*w(i,j,k,kbz)*time_now**2
        bb0=g*vdotB
        bb1=w(i,j,k,kbx)/g+g*vdotB*w(i,j,k,kvx)
        bb2=w(i,j,k,kby)/g+g*vdotB*w(i,j,k,kvy)
        bb3=w(i,j,k,kbz)/g+g*vdotB*w(i,j,k,kvz)
        !eq A6 of https://arxiv.org/pdf/1904.01530.pdf
        v(i,j,k,krc)=-(-om0*bb0+om1*bb1+om2*bb2+om3*bb3*time_now**2)
       end do   
    end do
  end do
 end if

else !2D case

 if(maxwell) then
  do k=iz1,iz2
    do j=iy1,iy2
      do i=ix1,ix2
        dBxcov_dy=(w(i,j+1,k,kbx)-w(i,j-1,k,kbx))/dy
        dBxcov_dz=0.
        dBycov_dx=(w(i+1,j,k,kby)-w(i-1,j,k,kby))/dx
        dBycov_dz=0.
        dBzcov_dx=g_cov(3)*(w(i+1,j,k,kbz)-w(i-1,j,k,kbz))/dx
        dBzcov_dy=g_cov(3)*(w(i,j+1,k,kbz)-w(i,j-1,k,kbz))/dy
        Ex_now=-time_now*(v(i,j,k,kvy)*v(i,j,k,kbz)-v(i,j,k,kvz)*v(i,j,k,kby))
        Ex_now_p=-time_now*(w(i+1,j,k,kvy)*w(i+1,j,k,kbz)-w(i+1,j,k,kvz)*w(i+1,j,k,kby))
        Ex_now_m=-time_now*(w(i-1,j,k,kvy)*w(i-1,j,k,kbz)-w(i-1,j,k,kvz)*w(i-1,j,k,kby))
        Ey_now=-time_now*(v(i,j,k,kvz)*v(i,j,k,kbx)-v(i,j,k,kvx)*v(i,j,k,kbz))
        Ey_now_p=-time_now*(w(i,j+1,k,kvz)*w(i,j+1,k,kbx)-w(i,j+1,k,kvx)*w(i,j+1,k,kbz))
        Ey_now_m=-time_now*(w(i,j-1,k,kvz)*w(i,j-1,k,kbx)-w(i,j-1,k,kvx)*w(i,j-1,k,kbz))
        Ez_now=-(v(i,j,k,kvx)*v(i,j,k,kby)-v(i,j,k,kvy)*v(i,j,k,kbx))/time_now
        Ez_now_p=-(w(i,j,k+1,kvx)*w(i,j,k+1,kby)-w(i,j,k+1,kvy)*w(i,j,k+1,kbx))/time_now
        Ez_now_m=-(w(i,j,k-1,kvx)*w(i,j,k-1,kby)-w(i,j,k-1,kvy)*w(i,j,k-1,kbx))/time_now
        Ex_past=-time_past*(vold(i,j,k,kvy)*vold(i,j,k,kbz)-vold(i,j,k,kvz)*vold(i,j,k,kby))
        Ey_past=-time_past*(vold(i,j,k,kvz)*vold(i,j,k,kbx)-vold(i,j,k,kvx)*vold(i,j,k,kbz))
        Ez_past=-(vold(i,j,k,kvx)*vold(i,j,k,kby)-vold(i,j,k,kvy)*vold(i,j,k,kbx))/time_past
        dExcont_dx=(Ex_now_p-Ex_now_m)/dx
        dEycont_dy=(Ey_now_p-Ey_now_m)/dy
        dEzcont_dz=0.
        dExcont_dt=(time_now*Ex_now-time_past*Ex_past)/dt
        dEycont_dt=(time_now*Ey_now-time_past*Ey_past)/dt
        dEzcont_dt=(time_now*Ez_now-time_past*Ez_past)/dt
        g=1/sqrt(1-v(i,j,k,kvx)**2-v(i,j,k,kvy)**2-v(i,j,k,kvz)**2*g_cov(3))
        divE=dExcont_dx+dEycont_dy+dEzcont_dz
        curlBx=dBzcov_dy-dBycov_dz
        curlBy=dBxcov_dz-dBzcov_dx
        curlBz=dBycov_dx-dBxcov_dy
        Jx=(curlBx-dExcont_dt)/time_now
        Jy=(curlBy-dEycont_dt)/time_now
        Jz=(curlBz-dEzcont_dt)/time_now
        Jdotv=Jx*v(i,j,k,kvx)+Jy*v(i,j,k,kvy)+Jz*v(i,j,k,kvz)*g_cov(3)
        v(i,j,k,krc)=g*(divE-Jdotv)
       end do   
    end do
  end do
 else !we compute the charge density using the vorticity
  do k=iz1,iz2
    do j=iy1,iy2
      do i=ix1,ix2
        g=1.d0/sqrt(1.d0-w(i,j,k,kvx)**2-w(i,j,k,kvy)**2-w(i,j,k,kvz)**2*time_now**2)
        uu0=-g!u0 covariant
        uu1=g*w(i,j,k,kvx)
        uu2=g*w(i,j,k,kvy)
        uu3=g*w(i,j,k,kvz)*time_now**2
        g_past=1.d0/sqrt(1.d0-vold(i,j,k,kvx)**2-vold(i,j,k,kvy)**2-vold(i,j,k,kvz)**2*time_past**2)
        uu0_past=-g_past!u0 covariant at previous timestep
        uu1_past=g_past*vold(i,j,k,kvx)
        uu2_past=g_past*vold(i,j,k,kvy)
        uu3_past=g_past*vold(i,j,k,kvz)*time_past**2
        grx=1.d0/sqrt(1.d0-w(i+1,j,k,kvx)**2-w(i+1,j,k,kvy)**2-w(i+1,j,k,kvz)**2*time_now**2)
        glx=1.d0/sqrt(1.d0-w(i-1,j,k,kvx)**2-w(i-1,j,k,kvy)**2-w(i-1,j,k,kvz)**2*time_now**2)
        gry=1.d0/sqrt(1.d0-w(i,j+1,k,kvx)**2-w(i,j+1,k,kvy)**2-w(i,j+1,k,kvz)**2*time_now**2)
        gly=1.d0/sqrt(1.d0-w(i,j-1,k,kvx)**2-w(i,j-1,k,kvy)**2-w(i,j-1,k,kvz)**2*time_now**2)
        u1ry=gry*w(i,j+1,k,kvx)
        u1ly=gly*w(i,j-1,k,kvx)
        u2rx=grx*w(i+1,j,k,kvy)
        u2lx=glx*w(i-1,j,k,kvy)
        u3rx=grx*w(i+1,j,k,kvz)
        u3lx=glx*w(i-1,j,k,kvz)
        u3ry=gry*w(i,j+1,k,kvz)
        u3ly=gly*w(i,j-1,k,kvz)
        !we recall that dx, dy and dz are already two times the cell spacing, see their definiton above
        du10=(uu1-uu1_past)/dt!du_x/dt
        du20=(uu2-uu2_past)/dt!du_y/dt
        du30=(uu3-uu3_past)/dt!du_z/dt
        du01=-(grx-glx)/dx!du_0/dx=-d_gammaLorentz/dx
        du02=-(gry-gly)/dy!du_0/dy=-d_gammaLorentz/dy
        du03=0.d0!du_0/dz
        du12=(u1ry-u1ly)/dy!du_x/dy
        du13=0.d0!du_x/dz
        du21=(u2rx-u2lx)/dx!du_y/dx
        du23=0.d0!du_y/dz
        du31=(u3rx-u3lx)/dx!du_z/dx
        du32=(u3ry-u3ly)/dy!du_z/dy
        !we use eq. A1 in https://arxiv.org/pdf/1904.01530.pdf
        om0=(e1203*du21*uu3+e2103*du12*uu3+e3102*du13*uu2+e1302*du31*uu2+e2301*du32*uu1+e3201*du23*uu1)/time_now
        om1=(e0213*du20*uu3+e2013*du02*uu3+e0312*du30*uu2+e3012*du03*uu2+e2310*du32*uu0+e3210*du23*uu0)/time_now
        om2=(e0123*du10*uu3+e1023*du01*uu3+e0321*du30*uu1+e3021*du03*uu1+e1320*du31*uu0+e3120*du13*uu0)/time_now
        om3=(e0231*du20*uu1+e2031*du02*uu1+e0132*du10*uu2+e1032*du01*uu2+e1230*du21*uu0+e2130*du12*uu0)/time_now
        vdotB=w(i,j,k,kvx)*w(i,j,k,kbx)+w(i,j,k,kvy)*w(i,j,k,kby)+w(i,j,k,kvz)*w(i,j,k,kbz)*time_now**2
        bb0=g*vdotB
        bb1=w(i,j,k,kbx)/g+g*vdotB*w(i,j,k,kvx)
        bb2=w(i,j,k,kby)/g+g*vdotB*w(i,j,k,kvy)
        bb3=w(i,j,k,kbz)/g+g*vdotB*w(i,j,k,kvz)
        !eq A6 of https://arxiv.org/pdf/1904.01530.pdf
        v(i,j,k,krc)=-(-om0*bb0+om1*bb1+om2*bb2+om3*bb3*time_now**2)
       end do   
    end do
  end do
 end if

end if

time_comp_rho_comov=t

end subroutine compute_rho_comoving_frame

! *****************************************************************************

subroutine out_derivatives(iout, compute_or_print)
!-- Write data on file 'der000n.dat'
  use system_eqgp, only: g_cov
  use work
  use evolve
  use viscous
  
  integer,intent(inout) :: iout
  integer      :: lv,i,lvm,lvp,proc
  integer      :: iii,jjj,kkk,lll,mmm,nnn
  character(len=7) :: fname
  integer :: ix,iy,iz

  integer :: error_flag, allocate_result
  
  
  real(8) :: v2, glf, vx, vy, vz, sum_der, sum_der_t, sum_der_0, sum_der_x, sum_der_y, sum_der_z
  real(8) :: dutdt, dutdx, dutdy, dutdz, duxdx, duydy, duzdz, duzdt, duzdx, duzdy, duxdt, duxdy, duxdz, duydt, duydx, duydz
  real(8) :: dtedt, dtedx, dtedy, dtedz, temperature

  real(8) :: uxr, uyr, uzr, uxl, uyl, uzl, glfr, glfl, tempr, templ, ds, enr, enl

  integer :: ixc, iyc, izc !used to simplify some indexese inside a loop

  !bitbucket variables just to call correctly the get_derived_data subroutine
  real(8) :: energy_density_bb, entropy_density_bb, temperature_bb

  integer, parameter :: kind_of_deriv=0 !if 1 it uses the reconstruction algorithm, otherwise the second order derivative 

  integer, intent(in) :: compute_or_print !if 0 it only computes derivatives, if 1 it also prints them

  error_flag=0

  if (pe0 .and. (compute_or_print .eq. 1)) then
    if (iout<10000) write (fname,'(a3,i4)') 'der' ,iout
    if (iout<1000) write (fname,'(a4,i3)') 'der0' ,iout
    if (iout<100) write (fname,'(a5,i2)') 'der00' ,iout
    if (iout<10) write (fname,'(a6,i1)') 'der000', iout

    if (id_of_the_run<10000) write (dir_out,'(a4,i4,a1)') 'outr' ,id_of_the_run,'/'
    if (id_of_the_run<1000) write (dir_out,'(a5,i3,a1)') 'outr0' ,id_of_the_run,'/'
    if (id_of_the_run<100) write (dir_out,'(a6,i2,a1)') 'outr00' ,id_of_the_run,'/'
    if (id_of_the_run<10) write (dir_out,'(a7,i1,a1)') 'outr000', id_of_the_run,'/'
    if(commandline_outdir_length .gt. 0) then
      prefix_dir=trim(adjustl(custom_output_directory))//"/"//dir_out
    else
      prefix_dir=dir_out
    end if
       
    open  (15,file=prefix_dir//fname//'.dat',form='unformatted', access='stream')
    if(output_precision .eq. 8) then
       write (15) t
    else
       write (15) real(t,4)
    end if

  end if
  
  w(ix1:ix2,iy1:iy2,iz1:iz2,1:nv)=v(ix1:ix2,iy1:iy2,iz1:iz2,1:nv)

  if(kind_of_deriv .eq. 1) then
    if (nx>1) then
      call evolve_bcx(ix1,iy1,iz1,1,nv,kv(1:nv))
    
      do iz=iz1,iz2
       do iy=iy1,iy2
        do ix=1-ngc, mx+ngc
        
          ixc=ix1+ix-1
     
          wd(ix,1)=1./sqrt(1.-w(ixc,iy,iz,kvx)**2.-w(ixc,iy,iz,kvy)**2.-g_cov(3)*w(ixc,iy,iz,kvz)**2.)
          wd(ix,2)=wd(ix,1)*w(ixc,iy,iz,kvx)
          wd(ix,3)=wd(ix,1)*w(ixc,iy,iz,kvy)
          wd(ix,4)=wd(ix,1)*w(ixc,iy,iz,kvz)
          if(freeze_type .eq. 0) then
            call get_derived_data(w(ixc,iy,iz,krh), w(ixc,iy,iz,kpr), energy_density_bb, wd(ix,5), entropy_density_bb,error_flag)
          else
            call get_derived_data(w(ixc,iy,iz,krh), w(ixc,iy,iz,kpr), wd(ix,5), temperature_bb, entropy_density_bb,error_flag)
          end if
          if(error_flag .gt. 0) then
            write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
            errcode=error_flag
            write(*,*) "Error code:", errcode
            write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
            run_crashed=.true.
            return
          end if
      
        end do
       
        call work_rec(mx,1,5,1)
      
        deriv(ix1:ix2,iy,iz,dtx)=ddx(ix1:ix2)*0.5*( wdl(1:mx,1) + wdr(1:mx,1) - wdl(0:mx-1,1) - wdr(0:mx-1,1))
        deriv(ix1:ix2,iy,iz,dxx)=ddx(ix1:ix2)*0.5*( wdl(1:mx,2) + wdr(1:mx,2) - wdl(0:mx-1,2) - wdr(0:mx-1,2))
        deriv(ix1:ix2,iy,iz,dyx)=ddx(ix1:ix2)*0.5*( wdl(1:mx,3) + wdr(1:mx,3) - wdl(0:mx-1,3) - wdr(0:mx-1,3))
        deriv(ix1:ix2,iy,iz,dzx)=ddx(ix1:ix2)*0.5*( wdl(1:mx,4) + wdr(1:mx,4) - wdl(0:mx-1,4) - wdr(0:mx-1,4))
        deriv(ix1:ix2,iy,iz,dtex)=ddx(ix1:ix2)*0.5*( wdl(1:mx,5) + wdr(1:mx,5) - wdl(0:mx-1,5) - wdr(0:mx-1,5))
       end do
      end do
    end if
   
    if(ny>1) then
      call evolve_bcy(ix1,iy1,iz1,1,nv,kv(1:nv))
      do iz=iz1,iz2
       do ix=ix1,ix2
        do iy=1-ngc, my+ngc
       
         iyc=iy1+iy-1

         wd(iy,1)=1./sqrt(1.-w(ix,iyc,iz,kvx)**2.-w(ix,iyc,iz,kvy)**2.-g_cov(3)*w(ix,iyc,iz,kvz)**2.)
         wd(iy,2)=wd(iy,1)*w(ix,iyc,iz,kvx)
         wd(iy,3)=wd(iy,1)*w(ix,iyc,iz,kvy)
         wd(iy,4)=wd(iy,1)*w(ix,iyc,iz,kvz)
         if(freeze_type .eq. 0) then
           call get_derived_data(w(ix,iyc,iz,krh), w(ix,iyc,iz,kpr), energy_density_bb, wd(iy,5), entropy_density_bb,error_flag)
         else
           call get_derived_data(w(ix,iyc,iz,krh), w(ix,iyc,iz,kpr), wd(iy,5), temperature_bb, entropy_density_bb,error_flag)
         end if
         if(error_flag .gt. 0) then
           write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
           errcode=error_flag
           write(*,*) "Error code:", errcode
           write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
           run_crashed=.true.
           return
         end if
        end do

        call work_rec(my,1,5,1)
      
        deriv(ix,iy1:iy2,iz,dty)=ddy(iy1:iy2)*0.5*( wdl(1:my,1) + wdr(1:my,1) - wdl(0:my-1,1) - wdr(0:my-1,1))
        deriv(ix,iy1:iy2,iz,dxy)=ddy(iy1:iy2)*0.5*( wdl(1:my,2) + wdr(1:my,2) - wdl(0:my-1,2) - wdr(0:my-1,2))
        deriv(ix,iy1:iy2,iz,dyy)=ddy(iy1:iy2)*0.5*( wdl(1:my,3) + wdr(1:my,3) - wdl(0:my-1,3) - wdr(0:my-1,3))
        deriv(ix,iy1:iy2,iz,dzy)=ddy(iy1:iy2)*0.5*( wdl(1:my,4) + wdr(1:my,4) - wdl(0:my-1,4) - wdr(0:my-1,4))
        deriv(ix,iy1:iy2,iz,dtey)=ddy(iy1:iy2)*0.5*( wdl(1:my,5) + wdr(1:my,5) - wdl(0:my-1,5) - wdr(0:my-1,5))
       end do
      end do
    end if
   
    if (nz>1) then

      call evolve_bcz(ix1,iy1,iz1,1,nv,kv(1:nv))
      do iy=iy1,iy2
       do ix=ix1,ix2
        do iz=1-ngc, mz+ngc
       
         izc=iz1+iz-1

         wd(iz,1)=1./sqrt(1.-w(ix,iy,izc,kvx)**2.-w(ix,iy,izc,kvy)**2.-g_cov(3)*w(ix,iy,izc,kvz)**2.)
         wd(iz,2)=wd(iz,1)*w(ix,iy,izc,kvx)
         wd(iz,3)=wd(iz,1)*w(ix,iy,izc,kvy)
         wd(iz,4)=wd(iz,1)*w(ix,iy,izc,kvz)
         if(freeze_type .eq. 0) then
           call get_derived_data(w(ix,iy,izc,krh), w(ix,iy,izc,kpr), energy_density_bb, wd(iz,5), entropy_density_bb,error_flag)
         else
           call get_derived_data(w(ix,iy,izc,krh), w(ix,iy,izc,kpr), wd(iz,5), temperature_bb, entropy_density_bb,error_flag)
         end if
         if(error_flag .gt. 0) then
           write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
           errcode=error_flag
           write(*,*) "Error code:", errcode
           write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
           run_crashed=.true.
           return
         end if
      
        end do

        call work_rec(mz,1,5,1)
      
        deriv(ix,iy,iz1:iz2,dtz)=ddz(iz1:iz2)*0.5*( wdl(1:mz,1) + wdr(1:mz,1) - wdl(0:mz-1,1) - wdr(0:mz-1,1))
        deriv(ix,iy,iz1:iz2,dxz)=ddz(iz1:iz2)*0.5*( wdl(1:mz,2) + wdr(1:mz,2) - wdl(0:mz-1,2) - wdr(0:mz-1,2))
        deriv(ix,iy,iz1:iz2,dyz)=ddz(iz1:iz2)*0.5*( wdl(1:mz,3) + wdr(1:mz,3) - wdl(0:mz-1,3) - wdr(0:mz-1,3))
        deriv(ix,iy,iz1:iz2,dzz)=ddz(iz1:iz2)*0.5*( wdl(1:mz,4) + wdr(1:mz,4) - wdl(0:mz-1,4) - wdr(0:mz-1,4))
        deriv(ix,iy,iz1:iz2,dtez)=ddz(iz1:iz2)*0.5*( wdl(1:mz,5) + wdr(1:mz,5) - wdl(0:mz-1,5) - wdr(0:mz-1,5))
       end do
      end do
    end if
  else !kind of deriv is not 1, so we compute derivatives with second order accuracy procedure
    if (nx>1) then
      call evolve_bcx(ix1,iy1,iz1,1,nv,kv(1:nv))
    
      do iz=iz1,iz2
       do iy=iy1,iy2
        do ix=ix1,ix2

          ds=1.d0/ddx(ix)+0.5d0/ddx(ix+1)+0.5d0/ddx(ix-1)
        
          glfr=1./sqrt(1.-w(ix+1,iy,iz,kvx)**2.-w(ix+1,iy,iz,kvy)**2.-g_cov(3)*w(ix+1,iy,iz,kvz)**2.)
          glfl=1./sqrt(1.-w(ix-1,iy,iz,kvx)**2.-w(ix-1,iy,iz,kvy)**2.-g_cov(3)*w(ix-1,iy,iz,kvz)**2.)
          uxr=glfr*w(ix+1,iy,iz,kvx)
          uyr=glfr*w(ix+1,iy,iz,kvy)
          uzr=glfr*w(ix+1,iy,iz,kvz)
          uxl=glfl*w(ix-1,iy,iz,kvx)
          uyl=glfl*w(ix-1,iy,iz,kvy)
          uzl=glfl*w(ix-1,iy,iz,kvz)
          call get_derived_data(w(ix-1,iy,iz,krh), w(ix-1,iy,iz,kpr), enl, templ, entropy_density_bb,error_flag)
          call get_derived_data(w(ix+1,iy,iz,krh), w(ix+1,iy,iz,kpr), enr, tempr, entropy_density_bb,error_flag)
          if(error_flag .gt. 0) then
            write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
            errcode=error_flag
            write(*,*) "Error code:", errcode
            write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
            run_crashed=.true.
            return
          end if
      
       
          deriv(ix,iy,iz,dtx)=(glfr-glfl)/ds
          deriv(ix,iy,iz,dxx)=(uxr-uxl)/ds
          deriv(ix,iy,iz,dyx)=(uyr-uyl)/ds
          deriv(ix,iy,iz,dzx)=(uzr-uzl)/ds
          if(freeze_type .eq. 0) then
            deriv(ix,iy,iz,dtex)=(tempr-templ)/ds
          else
            deriv(ix,iy,iz,dtex)=(enr-enl)/ds
          end if
        end do
       end do
      end do
    end if

    if(ny>1) then
      call evolve_bcy(ix1,iy1,iz1,1,nv,kv(1:nv))
      do iz=iz1,iz2
       do ix=ix1,ix2
        do iy=iy1,iy2

          ds=1./ddy(iy)+0.5/ddy(iy+1)+0.5/ddy(iy-1)
        
          glfr=1./sqrt(1.-w(ix,iy+1,iz,kvx)**2.-w(ix,iy+1,iz,kvy)**2.-g_cov(3)*w(ix,iy+1,iz,kvz)**2.)
          glfl=1./sqrt(1.-w(ix,iy-1,iz,kvx)**2.-w(ix,iy-1,iz,kvy)**2.-g_cov(3)*w(ix,iy-1,iz,kvz)**2.)
          uxr=glfr*w(ix,iy+1,iz,kvx)
          uyr=glfr*w(ix,iy+1,iz,kvy)
          uzr=glfr*w(ix,iy+1,iz,kvz)
          uxl=glfl*w(ix,iy-1,iz,kvx)
          uyl=glfl*w(ix,iy-1,iz,kvy)
          uzl=glfl*w(ix,iy-1,iz,kvz)
          call get_derived_data(w(ix,iy-1,iz,krh), w(ix,iy-1,iz,kpr), energy_density_bb, templ, entropy_density_bb,error_flag)
          call get_derived_data(w(ix,iy+1,iz,krh), w(ix,iy+1,iz,kpr), energy_density_bb, tempr, entropy_density_bb,error_flag)
          if(error_flag .gt. 0) then
            write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
            errcode=error_flag
            write(*,*) "Error code:", errcode
            write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
            run_crashed=.true.
            return
          end if
      
       
          deriv(ix,iy,iz,dty)=(glfr-glfl)/ds
          deriv(ix,iy,iz,dxy)=(uxr-uxl)/ds
          deriv(ix,iy,iz,dyy)=(uyr-uyl)/ds
          deriv(ix,iy,iz,dzy)=(uzr-uzl)/ds
          if(freeze_type .eq. 0) then
            deriv(ix,iy,iz,dtey)=(tempr-templ)/ds
          else
            deriv(ix,iy,iz,dtey)=(enr-enl)/ds
          end if

        end do
       end do
      end do
    end if
   
    if (nz>1) then

      call evolve_bcz(ix1,iy1,iz1,1,nv,kv(1:nv))
      do iy=iy1,iy2
       do ix=ix1,ix2
        do iz=iz1,iz2
       

          ds=1./ddz(iz)+0.5/ddz(iz+1)+0.5/ddz(iz-1)
        
          glfr=1./sqrt(1.-w(ix,iy,iz+1,kvx)**2.-w(ix,iy,iz+1,kvy)**2.-g_cov(3)*w(ix,iy,iz+1,kvz)**2.)
          glfl=1./sqrt(1.-w(ix,iy,iz-1,kvx)**2.-w(ix,iy,iz-1,kvy)**2.-g_cov(3)*w(ix,iy,iz-1,kvz)**2.)
          uxr=glfr*w(ix,iy,iz+1,kvx)
          uyr=glfr*w(ix,iy,iz+1,kvy)
          uzr=glfr*w(ix,iy,iz+1,kvz)
          uxl=glfl*w(ix,iy,iz-1,kvx)
          uyl=glfl*w(ix,iy,iz-1,kvy)
          uzl=glfl*w(ix,iy,iz-1,kvz)
          call get_derived_data(w(ix,iy,iz-1,krh), w(ix,iy,iz-1,kpr), energy_density_bb, templ, entropy_density_bb,error_flag)
          call get_derived_data(w(ix,iy,iz+1,krh), w(ix,iy,iz+1,kpr), energy_density_bb, tempr, entropy_density_bb,error_flag)
          if(error_flag .gt. 0) then
            write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
            errcode=error_flag
            write(*,*) "Error code:", errcode
            write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
            run_crashed=.true.
            return
          end if
      
       
          deriv(ix,iy,iz,dtz)=(glfr-glfl)/ds
          deriv(ix,iy,iz,dxz)=(uxr-uxl)/ds
          deriv(ix,iy,iz,dyz)=(uyr-uyl)/ds
          deriv(ix,iy,iz,dzz)=(uzr-uzl)/ds
          if(freeze_type .eq. 0) then
            deriv(ix,iy,iz,dtez)=(tempr-templ)/ds
          else
            deriv(ix,iy,iz,dtez)=(enr-enl)/ds
          end if
        end do
       end do
      end do
    end if
  end if !end if kind_of_deriv

  call dts(v,vold)
  
  do iz=iz1,iz2
    do iy=iy1,iy2
      do ix=ix1,ix2
         dutdt=deriv(ix,iy,iz,dtt)
         dutdx=deriv(ix,iy,iz,dtx)
         dutdy=deriv(ix,iy,iz,dty)
         dutdz=deriv(ix,iy,iz,dtz)
         duxdt=deriv(ix,iy,iz,dxt)
         duxdx=deriv(ix,iy,iz,dxx)
         duxdy=deriv(ix,iy,iz,dxy)
         duxdz=deriv(ix,iy,iz,dxz)
         duydt=deriv(ix,iy,iz,dyt)
         duydx=deriv(ix,iy,iz,dyx)
         duydy=deriv(ix,iy,iz,dyy)
         duydz=deriv(ix,iy,iz,dyz)
         duzdt=deriv(ix,iy,iz,dzt)
         duzdx=deriv(ix,iy,iz,dzx)
         duzdy=deriv(ix,iy,iz,dzy)
         duzdz=deriv(ix,iy,iz,dzz)
         dtedt=deriv(ix,iy,iz,dtet)
         dtedx=deriv(ix,iy,iz,dtex)
         dtedy=deriv(ix,iy,iz,dtey)
         dtedz=deriv(ix,iy,iz,dtez)
         vx=v(ix,iy,iz,kvx)
         vy=v(ix,iy,iz,kvy)
         vz=v(ix,iy,iz,kvz)

!DDD piece of code probably useless which can be safely removed
!         call get_derived_data(v(ix,iy,iz,krh),v(ix,iy,iz,kpr),energy_density_bb,temperature,entropy_density_bb,error_flag)
!         if(error_flag .gt. 0) then
!           write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
!           errcode=error_flag
!           write(*,*) "Error code:", errcode
!           write(*,*) "Position on the grid: ix=",lll,"iy=",mmm,"iz=",nnn
!           run_crashed=.true.
!           return
!         end if
         
      end do
    end do
  end do

    
  if(prl) then
    do i=1,nderivatives
     do nnn=1,mz
      do mmm=1,my
        call MPI_Gatherv(deriv(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,derivatives_all(1,mmm,nnn,i),recvcounts,displs,&
              &mpi_realtype,0,icomm,ierr)
      end do
     end do
    end do

  else
   derivatives_all=deriv
  end if
 
  if(pe0 .and. (compute_or_print .eq. 1)) then
   if(output_precision .eq. 8) then
    do nnn=1,nz
     do mmm=1,ny
      do lll=1,nx
       write(15) derivatives_all(lll,mmm,nnn,:)
      enddo
     enddo
    enddo
   else
    do nnn=1,nz
     do mmm=1,ny
      do lll=1,nx
        write(15) real(derivatives_all(lll,mmm,nnn,:),4)
      enddo
     enddo
    enddo
   end if
   close (15)
  end if
  
end subroutine out_derivatives

! *****************************************************************************

subroutine out_grid
!-- Write data on file 'grid.dat'

!  use parallel
!  use common
        implicit none
        integer qqq !just for counting
        logical dir_exists
        integer ierr
        
  if (pe0) then

    if (id_of_the_run<10000) write (dir_out,'(a4,i4,a1)') 'outr' ,id_of_the_run,'/'
    if (id_of_the_run<1000) write (dir_out,'(a5,i3,a1)') 'outr0' ,id_of_the_run,'/'
    if (id_of_the_run<100) write (dir_out,'(a6,i2,a1)') 'outr00' ,id_of_the_run,'/'
    if (id_of_the_run<10) write (dir_out,'(a7,i1,a1)') 'outr000', id_of_the_run,'/'

    if(commandline_outdir_length .gt. 0) then
      prefix_dir=trim(adjustl(custom_output_directory))//"/"//dir_out
    else
      prefix_dir=dir_out
    end if

    !check if the outpt directory can be created
    call EXECUTE_COMMAND_LINE('mkdir -p '//prefix_dir//'')
    write(*,*) 'mkdir '//dir_out//''
    call EXECUTE_COMMAND_LINE('touch '//prefix_dir//'testfile')
    inquire(file=prefix_dir//'/testfile',exist=dir_exists)
    if(dir_exists) then
      call EXECUTE_COMMAND_LINE('rm '//prefix_dir//'testfile')
    else
      write(*,*) 'Not able to create output directory, exiting...'
      write(*,*) 'Output directory should have been: ', prefix_dir
      call exit(1)
    end if

    open  (11,file=prefix_dir//'grid_summary.dat',status='REPLACE',form='formatted')

        write(11,*) nx,ny,nz
        write(11,*) d_ini
        write(11,*) xlim(1,:)
        write(11,*) xlim(2,:)

    close (11)

    open  (11,file=prefix_dir//'grid.dat',status='REPLACE',form='formatted')
    write (11,*) nx, ny, nz
    write (11,*) (x(qqq),qqq=1,nx)
    write (11,*) (y(qqq),qqq=1,ny)
    write (11,*) (z(qqq),qqq=1,nz)
    close (11)
  end if

  if(prl) call MPI_Barrier(MPI_COMM_WORLD,ierr)

end subroutine out_grid

! *****************************************************************************

subroutine out_dump(iout)
!-- Write data on file 'dump.dat'

  use hypersurface
  implicit none
  integer,intent(in) :: iout
  character(len=8) :: dumpfile
  integer :: i, file_status,mmm,nnn
  
  if(prl) then
    do i=1,nv
     do nnn=1,mz
      do mmm=1,my
        call MPI_Gatherv(u(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,uall(1,mmm,nnn,i),recvcounts,displs,mpi_realtype,0&
             &,icomm,ierr)
      end do
     end do
    end do
    do i=krh,kvold_end
     do nnn=1,mz
      do mmm=1,my
        call MPI_Gatherv(uold(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,uallold(1,mmm,nnn,i),recvcounts,displs,mpi_realtype,0&
             &,icomm,ierr)
        call MPI_Gatherv(vold(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,vallold(1,mmm,nnn,i),recvcounts,displs,mpi_realtype,0&
             &,icomm,ierr)
      end do
     end do
    end do
  else
    ! here we are copying the whole arrays, but we leave the final index to make clear what variables are inside the 4th index
    uall(:,:,:,1:nv)=u(:,:,:,1:nv)
    vallold(:,:,:,krh:kvold_end)=vold(:,:,:,krh:kvold_end)
    uallold(:,:,:,krh:kvold_end)=uold(:,:,:,krh:kvold_end)
  end if


  if (pe0) then

    if (id_of_the_run<10000) write (dir_out,'(a4,i4,a1)') 'outr' ,id_of_the_run,'/'
    if (id_of_the_run<1000) write (dir_out,'(a5,i3,a1)') 'outr0' ,id_of_the_run,'/'
    if (id_of_the_run<100) write (dir_out,'(a6,i2,a1)') 'outr00' ,id_of_the_run,'/'
    if (id_of_the_run<10) write (dir_out,'(a7,i1,a1)') 'outr000', id_of_the_run,'/'

    if(commandline_outdir_length .gt. 0) then
      prefix_dir=trim(adjustl(custom_output_directory))//"/"//dir_out
    else
      prefix_dir=dir_out
    end if

    open (10,file=prefix_dir//'dump.dat',form='unformatted',iostat=file_status, access='stream')
    if(file_status .ne. 0) then  
       write(*,*) "Sorry, but I can't open the file:"
       write(*,*) prefix_dir//'dump.dat'
       write(*,*) "So, I cannot write the dump.dat file for restarting (as chosen with option restart_type=1)"
       write(*,*) "I think it is better to quit..."
       call exit(1)
    end if
      

    !we didn't collect v_backup, vall, time_arr and der_backup because they are already on the 0 processor
    write (10) iout
    write (10) t
    write (10) timeold
    write (10) timeinterval
    write (10) vall
    write (10) uall
    write (10) vallold
    write (10) uallold
    write (10) g_cov_3_old
    if(out_freeze) then
      write (10) thyp
      write (10) time_arr
      write (10) v_backup
      write (10) field(2,:,:,:)
      if (derivatives_out) write (10) der_backup
    end if
    close (10)
  endif

end subroutine out_dump

! *****************************************************************************

subroutine out_restart(iout)
!-- It reads data fron file 'dump.dat'

  use hypersurface
  implicit none
  integer,intent(out) :: iout
  character(len=8) :: dumpfile
  integer :: i,file_status, mmm, nnn
        

  if (pe0) then
    if (id_of_the_run<10000) write (dir_out,'(a4,i4,a1)') 'outr' ,iout,'/'
    if (id_of_the_run<1000) write (dir_out,'(a5,i3,a1)') 'outr0' ,iout,'/'
    if (id_of_the_run<100) write (dir_out,'(a6,i2,a1)') 'outr00' ,iout,'/'
    if (id_of_the_run<10) write (dir_out,'(a7,i1,a1)') 'outr000', iout,'/'

    if(commandline_outdir_length .gt. 0) then
      prefix_dir=trim(adjustl(custom_output_directory))//"/"//dir_out
    else
      prefix_dir=dir_out
    end if

    open (10,file=prefix_dir//'dump.dat',form='unformatted',iostat=file_status, access='stream')
    if(file_status .ne. 0) then  
      write(*,*) "Sorry, but I can't open the file:"
      write(*,*) prefix_dir//'dump.dat'
      write(*,*) "So, I cannot restart from the last frame and I have to quit..."
      call exit(1)
    end if

    read (10) iout
    read (10) t
    read (10) timeold
    read (10) timeinterval
    read (10) vall
    read (10) uall
    read (10) vallold
    read (10) uallold
    read (10) g_cov_3_old
    if(out_freeze) then
      read (10) thyp
      read (10) time_arr
      read (10) v_backup
      read (10) field(2,:,:,:)
      if(derivatives_out) read (10) der_backup
    end if
    close (10)
  endif

  if (prl) then
    call mpi_bcast(iout,1,mpi_integer ,0,icomm,ierr)
    call mpi_bcast(t   ,1,mpi_realtype,0,icomm,ierr)
    call mpi_bcast(timeold,1,mpi_realtype,0,icomm,ierr)
    call mpi_bcast(timeinterval,1,mpi_realtype,0,icomm,ierr)
    call mpi_bcast(g_cov_3_old,1,mpi_realtype,0,icomm,ierr)
    do i=1,nv
     do nnn=1,mz
      do mmm=1,my
        call MPI_Scatterv(vall(1,mmm,nnn,i),recvcounts,displs,mpi_realtype,v(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,0&
             &,icomm,ierr)
        call MPI_Scatterv(uall(1,mmm,nnn,i),recvcounts,displs,mpi_realtype,u(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,0&
             &,icomm,ierr)
      end do
     end do
    end do
    do i=krh,kvold_end
     do nnn=1,mz
      do mmm=1,my
        call MPI_Scatterv(vallold(1,mmm,nnn,i),recvcounts,displs,mpi_realtype,vold(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,0&
             &,icomm,ierr)
        call MPI_Scatterv(uallold(1,mmm,nnn,i),recvcounts,displs,mpi_realtype,uold(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,0&
             &,icomm,ierr)
      end do
     end do
    end do

  call MPI_Barrier(icomm, ierr)
 
  else
    v=vall
    u=uall
    vold=vallold
    uold=uallold
  end if

  if(pe0 .and. derivatives_out) primitives=vall

end subroutine out_restart

! *****************************************************************************

subroutine check_temp()
  use eos
  implicit none
  integer :: nnn, mmm, lll
  real(8) :: energy_density_value, temperature_value, temp_max_prl, temp_max, en_max, en_max_prl
  integer :: errcode

  errcode=0 
  
  temp_max_prl=0.
  temp_max=0.
  en_max=0.
  en_max_prl=0.
  
  !reminder: freeze_type=0 means freeze out based on temperature, freeze_type=1 means freeze-out
  if(freeze_type .eq. 0) then
    do nnn=iz1,iz2
     do mmm=iy1,iy2
      do lll=ix1,ix2
       call get_derived_temp_for_check(v(lll,mmm,nnn,krh), v(lll,mmm,nnn,kpr), energy_density_value, temperature_value,errcode)
       temp_max=max(temp_max,temperature_value)
      end do
     end do
    end do

    if (prl) call MPI_Barrier(MPI_COMM_WORLD,ierr)

    if (prl) then
        call MPI_Reduce(temp_max, temp_max_prl, 1, mpi_realtype, MPI_MAX, 0, icomm, ierr)
    else
        temp_max_prl=temp_max
    end if

    if(pe0) then
      if(temp_max_prl .lt. temp_end) then
         temp_check=1
      else
         temp_check=0
      end if
    end if
  else
    do nnn=iz1,iz2
     do mmm=iy1,iy2
      do lll=ix1,ix2
       call get_derived_temp_for_check(v(lll,mmm,nnn,krh), v(lll,mmm,nnn,kpr), energy_density_value, temperature_value,errcode)
       en_max=max(en_max,energy_density_value)
      end do
     end do
    end do

    if (prl) call MPI_Barrier(MPI_COMM_WORLD,ierr)

    if (prl) then
        call MPI_Reduce(en_max, en_max_prl, 1, mpi_realtype, MPI_MAX, 0, icomm, ierr)
    else
        en_max_prl=en_max
    end if

    if(pe0) then
      !we remind that, when performing freezeout based on energy density, temp_end actually means en_end,
      !the energy density treshold to stop the simulation
      if(en_max_prl .lt. temp_end) then
         temp_check=1
      else
         temp_check=0
      end if
    end if
  end if
  if (prl) call mpi_bcast(temp_check,1,mpi_integer ,0,icomm,ierr)
  
end subroutine check_temp

! *****************************************************************************

subroutine out_copy_tools

  call EXECUTE_COMMAND_LINE('cp '//param_file//' '//prefix_dir)
  call EXECUTE_COMMAND_LINE('mv config_summary.dat '//prefix_dir)

end subroutine out_copy_tools
! *****************************************************************************

subroutine out_hyp(iout)
  use common, only: tout, thyp
  use hypersurface
  implicit none
  logical, save :: first_time
  integer ix, iy, iz
  real(8) temp,en,et
  integer errcode_local,error_flag,iout
  integer idx !just a counter

  !we use this trick to distinguish between first real runs and first call to this subroutine after restarting
  if(time_arr(2) .eq. 0) first_time=.true.

  if(print_rho_comov) call compute_rho_comoving_frame
  
  if(t_collect_d_a .ne. t) call collect_distributed_arrays

  if(derivatives_out) call out_derivatives(iout,0)

  if(pe0) then
 
    if(first_time) time_arr(2)=t

    if(freeze_type .eq. 0) then
      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
            call get_derived_data(vall(ix,iy,iz,krh),vall(ix,iy,iz,kpr), en,temp, et,errcode_local)
            if(errcode_local .gt. 0) then
              write(*,*) "Error into the out_hyp subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
              write(*,*) "Error code:", errcode_local
              write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
              run_crashed=.true.
              return
            end if
            field(1,ix,iy,iz)=field(2,ix,iy,iz)
            field(2,ix,iy,iz)=temp
          end do
        end do
      end do
    else if(freeze_type .eq. 1) then
      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
            call get_derived_data(vall(ix,iy,iz,krh),vall(ix,iy,iz,kpr), en,temp, et,errcode_local)
            if(errcode_local .gt. 0) then
              write(*,*) "Error into the out_hyp subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
              write(*,*) "Error code:", errcode_local
              write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
              run_crashed=.true.
              return
            end if
            field(1,ix,iy,iz)=field(2,ix,iy,iz)
            field(2,ix,iy,iz)=en
          end do
        end do
      end do
    else
      write(*,*) "Error, it is not clear if freezout is based on temperature or energy density..."
      write(*,*) 'Please, check your parameters (usually, param.dat) file...'
      call exit(1)
    end if
 
    if(first_time) then
      first_time=.false.
      primitives=vall
    else

      v_backup=primitives
      primitives=vall
      time_arr(1)=time_arr(2)
      time_arr(2)=t
      call find_hypersurface()
    end if

    if(derivatives_out) der_backup=derivatives_all

    
  end if

 
end subroutine out_hyp
! *****************************************************************************

subroutine compute_on_transverse_plane(iout,error_flag)
  !AAA here we are assuming that cells on the grid all have the same size
 
  use common
  implicit none
  real(8), allocatable, dimension(:),save::  xcm_num, ycm_num, cm_den, ecc_num, ecc_den, ell_num, ell_den, dfl_num, dfl_den
  real(8), allocatable, dimension(:),save ::  xcm_numvec, ycm_numvec, cm_denvec, ecc_numvec, ecc_denvec, ell_numvec, ell_denvec
  real(8), allocatable, dimension(:), save :: dfl_numvec, dfl_denvec
  integer :: i,j, ix, iy,iz, error_flag, iout
  logical, save :: first_time=.true.
  integer :: allocate_error
  real(8) :: dx, dy
  real(8) :: energy_l, energy_r, v0_l, v0_r
  character(len=10) :: epfile, ecfile, dffile

  if (iout<10000) then
     write (epfile,'(a2,i4,a4)') 'ep' ,iout,'.dat'
     write (ecfile,'(a2,i4,a4)') 'ec' ,iout,'.dat'
     write (dffile,'(a2,i4,a4)') 'df' ,iout,'.dat'
  end if
  if (iout<1000) then
     write (epfile,'(a3,i3,a4)') 'ep0' ,iout,'.dat'
     write (ecfile,'(a3,i3,a4)') 'ec0' ,iout,'.dat'
     write (dffile,'(a3,i3,a4)') 'df0' ,iout,'.dat'
  end if
  if (iout<100) then
     write (epfile,'(a4,i2,a4)') 'ep00' ,iout,'.dat'
     write (ecfile,'(a4,i2,a4)') 'ec00' ,iout,'.dat'
     write (dffile,'(a4,i2,a4)') 'df00' ,iout,'.dat'
  end if

  if (iout<10) then
     write (epfile,'(a5,i1,a4)') 'ep000' ,iout,'.dat'
     write (ecfile,'(a5,i1,a4)') 'ec000' ,iout,'.dat'
     write (dffile,'(a5,i1,a4)') 'df000' ,iout,'.dat'
  end if


  if(first_time) then
     first_time=.false.
     allocate(xcm_numvec(1:ny), ycm_numvec(1:ny), cm_denvec(1:ny), ecc_numvec(1:ny), ecc_denvec(1:ny), ell_numvec(1:ny),&
             &ell_denvec(1:ny),dfl_numvec(1:ny),dfl_denvec(1:ny),STAT=allocate_error)
     if(allocate_error .ne. 0) then
       write(*,*) "Error in compute_on_transeverse_plane subroutine (system.f08)"
       write(*,*) "It was not possible to allocate the first series of dynamic arrays with (1:ny) elements, so I quit..."
       call exit(1)
     end if
     allocate(xcm_num(1:nz), ycm_num(1:nz), cm_den(1:nz), ecc_num(1:nz), ecc_den(1:nz), ell_num(1:nz),&
             &ell_den(1:nz),dfl_num(1:nz),dfl_den(1:nz),STAT=allocate_error)
     if(allocate_error .ne. 0) then
       write(*,*) "Error in compute_on_transeverse_plane subroutine (system.f08)"
       write(*,*) "It was not possible to allocate the second series of dynamic arrays with (1:nz) elements, so I quit..."
       call exit(1)
     end if
  end if

  open  (20,file=prefix_dir//epfile,status='REPLACE',form='formatted')
  open  (21,file=prefix_dir//ecfile,status='REPLACE',form='formatted')
  open  (22,file=prefix_dir//dffile,status='REPLACE',form='formatted')
  
  error_flag=0
  xcm_num=0.
  ycm_num=0.
  cm_den=0.
  ecc_num=0.
  ecc_den=0.
  dfl_num=0.
  dfl_den=0.
  
  dx=(x2-x1)/nx
  dy=(y2-y1)/ny
  
  do iz=1,nz !big loop
    xcm_numvec=0.
    ycm_numvec=0.
    cm_denvec=0.
    ecc_numvec=0.
    ecc_denvec=0.
    ell_numvec=0.
    ell_denvec=0.
    dfl_numvec=0.
    dfl_denvec=0.

    do iy=1,ny
      call eos_energy(vall(1,iy,iz,krh),energy_r,vall(1,iy,iz,kpr),errcode)
      if(errcode/=0) then
        write(*,*) "Error when computing energy_r into compute_on_transverse_plane subroutine, inside system.f08"
        write(*,*) 'Call before x-cycle: ix,iy,iz:',1,iy,iz
        write(*,*) 'Error:',errcode
        run_crashed=.true.
        return
      end if

      v0_r = 1./sqrt(1.-vall(1,iy,iz,kvx)*vall(1,iy,iz,kvx)*g_cov(1)-vall(1,iy,iz,kvy)*vall(1,iy,iz,kvy)*g_cov(2)-&
          & vall(1,iy,iz,kvz)*vall(1,iy,iz,kvz)*g_cov(3))
  
      do ix=1,nx-1
        energy_l=energy_r
        v0_l=v0_r
        v0_r = 1./sqrt(1.-vall(ix+1,iy,iz,kvx)*vall(ix+1,iy,iz,kvx)*g_cov(1)-vall(ix+1,iy,iz,kvy)*vall(ix+1,iy,iz,kvy)*g_cov(2)-&
          & vall(ix+1,iy,iz,kvz)*vall(ix+1,iy,iz,kvz)*g_cov(3))
        call eos_energy(vall(ix+1,iy,iz,krh),energy_r,vall(ix+1,iy,iz,kpr),errcode)
        if(errcode/=0) then
          write(*,*) "Error when computing energy_r into compute_on_transverse_plane subroutine, inside system.f08"
          write(*,*) 'ix+1,iy,iz:',ix+1,iy,iz
          write(*,*) 'Error:',errcode
          run_crashed=.true.
          return
         end if
        !in these formulas we multiply for dx/6 even if this term, appearing both in numerator and denominator, will be canceled 
        xcm_numvec(iy)=xcm_numvec(iy)+dx*(x(ix)*energy_l+x(ix+1)*energy_r)/2.
        ecc_numvec(iy)=ecc_numvec(iy)+((y(iy)*y(iy)-x(ix)*x(ix))*energy_l+(y(iy)*y(iy)-x(ix+1)*x(ix+1))*energy_r)/2.
        ell_numvec(iy)=ell_numvec(iy)+((energy_l+vall(ix,iy,iz,kpr))*(v0_l*v0_l)*(vall(ix,iy,iz,kvx)**2.-vall(ix,iy,iz,kvy)**2.)&
                      &+(energy_r+vall(ix+1,iy,iz,kpr))*(v0_r*v0_r)*(vall(ix+1,iy,iz,kvx)**2.-vall(ix+1,iy,iz,kvy)**2.))/2.
        ell_denvec(iy)=ell_denvec(iy)+((energy_l+vall(ix,iy,iz,kpr))*(v0_l*v0_l)*(vall(ix,iy,iz,kvx)**2.+vall(ix,iy,iz,kvy)**2.+&
                      &2.*vall(ix,iy,iz,kpr))+(energy_r+vall(ix+1,iy,iz,kpr))*(v0_r*v0_r)*(vall(ix+1,iy,iz,kvx)**2.+&
                      &vall(ix+1,iy,iz,kvy)**2.+2.*vall(ix+1,iy,iz,kpr)))/2.
        ycm_numvec(iy)=ycm_numvec(iy)+dx*(energy_l+energy_r)/2.
        !since we multiply for y only when integrating along y, for now ynumvec and denvec have the same values
        cm_denvec(iy)=cm_denvec(iy)+dx*(energy_l+energy_r)/2.
        ecc_denvec(iy)=ecc_denvec(iy)+((y(iy)*y(iy)+x(ix)*x(ix))*energy_l+(y(iy)*y(iy)+x(ix+1)*x(ix+1))*energy_r)/2.
        dfl_numvec(iy)=dfl_numvec(iy)+(vall(ix,iy,iz,kvx)*v0_l*energy_l+vall(ix+1,iy,iz,kvx)*v0_r*energy_r)/2.
        dfl_denvec(iy)=dfl_denvec(iy)+(v0_l*energy_l+v0_r*energy_r)/2.
      end do
    end do
 
    do iy=1,ny-1
      xcm_num(iz)=xcm_num(iz)+dy*(xcm_numvec(iy)+xcm_numvec(iy+1))/2.
      ycm_num(iz)=ycm_num(iz)+dy*(y(iy)*ycm_numvec(iy)+y(iy+1)*ycm_numvec(iy+1))/2.
      cm_den(iz)=cm_den(iz)+dy*(cm_denvec(iy)+cm_denvec(iy+1))/2.
      ecc_num(iz)=ecc_num(iz)+dy*(ecc_numvec(iy)+ecc_numvec(iy+1))/2.
      ecc_den(iz)=ecc_den(iz)+dy*(ecc_denvec(iy)+ecc_denvec(iy+1))/2.
      ell_num(iz)=ell_num(iz)+dy*(ell_numvec(iy)+ell_numvec(iy+1))/2.
      ell_den(iz)=ell_den(iz)+dy*(ell_denvec(iy)+ell_denvec(iy+1))/2.
      dfl_num(iz)=dfl_num(iz)+dy*(dfl_numvec(iy)+dfl_numvec(iy+1))/2.
      dfl_den(iz)=dfl_den(iz)+dy*(dfl_denvec(iy)+dfl_denvec(iy+1))/2.
      
    end do

    if(cm_den(iz) .eq. 0) then
      write(*,*) "Error in subroutine compute_on_transevers_plane (file system.f08)"
      write(*,*) "Variable cm_den=0 for iz=",iz
      write(*,*) "Maybe energy density is 0 all over the grid. Leaving..."
      error_flag=222
      run_crashed=.true.
      return
    end if
 
    if(ecc_den(iz) .eq. 0) then
      write(*,*) "Error in subroutine find_cm_on_transevers_plane (file system.f08)"
      write(*,*) "Variable ecc_den(iz)=0 for iz=",iz
      write(*,*) "Maybe density is 0 all over the grid. Leaving..."
      error_flag=223
      run_crashed=.true.
      return
    end if

    if(ell_den(iz) .eq. 0) then
      write(*,*) "Error in subroutine find_cm_on_transevers_plane (file system.f08)"
      write(*,*) "Variable ell_den(iz)=0 for iz=",iz
      write(*,*) "Maybe density is 0 all over the grid. Leaving..."
      error_flag=224
      run_crashed=.true.
      return
    end if

    if(dfl_den(iz) .eq. 0) then
      write(*,*) "Error in subroutine find_cm_on_transevers_plane (file system.f08)"
      write(*,*) "Variable dfl_den(iz)=0 for iz=",iz
      write(*,*) "Maybe density is 0 all over the grid. Leaving..."
      error_flag=225
      run_crashed=.true.
      return
    end if
 
 
    xcm(iz)=xcm_num(iz)/cm_den(iz)
    ycm(iz)=ycm_num(iz)/cm_den(iz)
    eccentricity(iz)=ecc_num(iz)/ecc_den(iz)
    elliptic_flow(iz)=ell_num(iz)/ell_den(iz)
    directed_flow(iz)=dfl_num(iz)/dfl_den(iz)
   
    write(20,*) z(iz), elliptic_flow(iz)
    write(21,*) z(iz), eccentricity(iz)
    write(22,*) z(iz), directed_flow(iz)
  

 end do ! end big loop on z
 close(20)
 close(21)
 close(22)

end subroutine compute_on_transverse_plane

! *************************************************************

subroutine collect_distributed_arrays
  integer      :: i,lll,mmm,nnn

  t_collect_d_a=t

  if(prl) then
    do i=1,nv
     do nnn=1,mz
      do mmm=1,my
        call MPI_Gatherv(v(ix1,iy1+mmm-1,iz1+nnn-1,i),mx,mpi_realtype,vall(1,mmm,nnn,i),recvcounts,displs,mpi_realtype,0&
             &,icomm,ierr)
      end do
     end do
    end do

    if(viscosity) then
     do nnn=1,mz
      do mmm=1,my
        call MPI_Gatherv(deriv(ix1,iy1+mmm-1,iz1+nnn-1,dtt),mx,mpi_realtype,vderivatives(1,mmm,nnn,1),recvcounts,displs,&
             &mpi_realtype,0,icomm,ierr)
        call MPI_Gatherv(deriv(ix1,iy1+mmm-1,iz1+nnn-1,dxx),mx,mpi_realtype,vderivatives(1,mmm,nnn,2),recvcounts,displs,&
             &mpi_realtype,0,icomm,ierr)
        call MPI_Gatherv(deriv(ix1,iy1+mmm-1,iz1+nnn-1,dyy),mx,mpi_realtype,vderivatives(1,mmm,nnn,3),recvcounts,displs,&
             &mpi_realtype,0,icomm,ierr)
        call MPI_Gatherv(deriv(ix1,iy1+mmm-1,iz1+nnn-1,dzz),mx,mpi_realtype,vderivatives(1,mmm,nnn,4),recvcounts,displs,&
             &mpi_realtype,0,icomm,ierr)
      end do
     end do
    end if
  else
   vall=v
   if(viscosity) then
     vderivatives(:,:,:,1)=deriv(:,:,:,dtt)
     vderivatives(:,:,:,2)=deriv(:,:,:,dxx)
     vderivatives(:,:,:,3)=deriv(:,:,:,dyy)
     vderivatives(:,:,:,4)=deriv(:,:,:,dzz)
    end if
  end if

end subroutine collect_distributed_arrays

! *************************************************************

end module out

