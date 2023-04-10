! *****************************************************************************
! *                                                                           * 
! *  ECHO-QGP                                                                 * 
! *                                                                           * 
! *  Version: 1.0.00                                                          *
! *                                                                           *
! *  Date: May, 22 - 2015                                                     *
! *                                                                           *
! *  File: readx.f90                                                          *
! *                                                                           *
! *  License: GPL version 2.0 (Please, read the file LICENSE.TXT)             *
! *                                                                           *
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
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
! *  Authors: Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)            *
! *                                                                           *
! *  Contributors: Andrea Beraudo (beraudo@to.infn.it)                        *
! *                                                                           *
! *  Acknowledgments:                                                         *
! *                                                                           *
! *****************************************************************************

      program readx
      implicit none
      
      integer neta, nx, ny
      integer ieta, ix, iy
      integer filerror
      integer ntstart,ntend,it,nv0,nv1,nv2,nv3
      integer out_precision
      integer i
      real etastep, xstep, ystep
      real vz
      real time
      real taustep
      real, allocatable, dimension(:,:,:) :: rhomat, vxmat, vymat, vetamat,prmat, epsmat, entropymat,tempmat,bulkmat
      real, allocatable, dimension(:,:,:) :: ttmat, txmat, tymat, tzmat, xymat, xzmat, yzmat, xxmat, yymat, zzmat, v0mat
      real, allocatable, dimension(:,:,:) :: dutdtmat, duxdxmat, duydymat, duzdzmat, thetamat
      real, allocatable, dimension(:) :: xvector, yvector, etavector
      real(KIND=4) etastep4, xstep4, ystep4
      real(KIND=4) vz4
      real(KIND=4) time4
      real(KIND=4) taustep4,time04
      real(KIND=4), allocatable, dimension(:,:,:) :: rhomat4, vxmat4, vymat4, vetamat4,prmat4, epsmat4, entropymat4,tempmat4
      real(KIND=4), allocatable, dimension(:,:,:) :: bulkmat4, ttmat4, txmat4, tymat4, tzmat4, xymat4, xzmat4, yzmat4, xxmat4 
      real(KIND=4), allocatable, dimension(:,:,:) :: yymat4, zzmat4, v0mat4,dutdtmat4,duxdxmat4,duydymat4, duzdzmat4,thetamat4
      character(len=4) ntplotname,no
      character(len=40) tempfile, vxfile, vyfile, vzfile, epsfile, rhofile, prfile, entropyfile,bulkfile
      character(len=40) ttfile, txfile, tyfile, tzfile, xyfile, xzfile, yzfile, xxfile, yyfile, zzfile, v0file
      character(len=40) dutdtfile, duxdxfile, duydyfile, duzdzfile, thetafile
      character(len=255) outdir
      LOGICAL :: directory_exists !to inquire if the output directory exists
      type out_type
       integer density
       integer vx
       integer vy
       integer vz
       integer pressure
       integer energy_density
       integer temperature
       integer entropy_density
       integer bulk
       integer pitt
       integer pitx
       integer pity
       integer pitz
       integer pixy
       integer pixz
       integer piyz
       integer pixx
       integer piyy
       integer pizz
       integer v0
       integer dutdt
       integer duxdx
       integer duydy
       integer duzdz
       integer theta
      end type out_type

      type(out_type) out_sel
      
 008 format(i1,35x)
 009 format(10x,f14.8)
      if (command_argument_count() .lt. 2) then
      write(*,*) 'Please, insert the range of echo-qgp output files from which you want to extract values.'
      write(*,*) '(Optionally, you can also write the name of the directory where the output files will be written)'
      write(*,*) '(default: postproc/readx)'
      call exit(1)
      end if
          
      !we get the starting index 
      call get_command_argument(1,ntplotname)
      !we convert from char to integer and we put the value in the variable ntstart
      read(ntplotname,*) ntstart
      !we get the starting index 
      call get_command_argument(2,ntplotname)
      !we convert from char to integer and we put the value in the variable ntstart
      read(ntplotname,*) ntend
      ! if we've got two arguments, we use the second one as the name of the directory where to write output files
      if (command_argument_count() .ge. 3) then ! the arguments after the second one are ignored
      call get_command_argument(3,outdir)
      else
      outdir='postproc/readx'
      endif
      
      ! check if the output directory exists      
      INQUIRE(FILE=trim(adjustl(outdir)), EXIST=directory_exists)
        if ( .NOT. directory_exists ) then
           call system('mkdir -p '//trim(adjustl(outdir)))
        endif
       
      open (2, file="grid_summary.dat", status='OLD', iostat=filerror)
      if (filerror .ne. 0) then
         write(*,*) "Sorry, but I need the file grid_summary.dat, which I cannot open, so I'm forced to quit!"
         close(2)
         call exit(2)
      end if

      read(2,*) nx, ny, neta
      read(2,*) xstep, ystep, etastep
      close(2)

      open (2, file="param.dat", status='OLD', iostat=filerror)
      if (filerror .ne. 0) then
         write(*,*) "Sorry, but I need the file param.dat, which I cannot open, so I'm forced to quit!"
         close(2)
         call exit(2)
      end if

      do i=1,24
         read(2,*) 
      end do
      read(2,009) taustep
      close(2)

      open (11, file="vars_summary.dat", status='OLD', iostat=filerror)
      if (filerror .ne. 0) then
         write(*,*) "Sorry, but I need the file vars_summary.dat, which I cannot open, so I'm forced to quit!"
         close(11)
         call exit(2)
      end if
     

      ! now we read the vars_summary.dat file to know which variables are contained in the output
      read(11,*)
      read(11,008) out_precision
      read(11,008) out_sel%density
      read(11,008) out_sel%vx
      read(11,008) out_sel%vy
      read(11,008) out_sel%vz
      read(11,008) out_sel%pressure
      read(11,008) out_sel%energy_density
      read(11,008) out_sel%temperature
      read(11,008) out_sel%entropy_density
      read(11,008) out_sel%bulk
      read(11,008) out_sel%pitt
      read(11,008) out_sel%pitx
      read(11,008) out_sel%pity
      read(11,008) out_sel%pitz
      read(11,008) out_sel%pixy
      read(11,008) out_sel%pixz
      read(11,008) out_sel%piyz
      read(11,008) out_sel%pixx
      read(11,008) out_sel%piyy
      read(11,008) out_sel%pizz
      read(11,008) out_sel%v0
      read(11,008) out_sel%dutdt
      read(11,008) out_sel%duxdx
      read(11,008) out_sel%duydy
      read(11,008) out_sel%duzdz
      read(11,008) out_sel%theta

      close(11)
     
      if(out_precision .eq. 8) then
        if(out_sel%density .eq. 1) allocate(rhomat(nx,ny,neta))
        if(out_sel%vx .eq. 1) allocate(vxmat(nx,ny,neta))
        if(out_sel%vy .eq. 1) allocate(vymat(nx,ny,neta))
        if(out_sel%vz .eq. 1) allocate(vetamat(nx,ny,neta))
        if(out_sel%pressure .eq. 1) allocate(prmat(nx,ny,neta))
        if(out_sel%energy_density .eq. 1) allocate(epsmat(nx,ny,neta))
        if(out_sel%temperature .eq. 1) allocate(tempmat(nx,ny,neta))
        if(out_sel%entropy_density .eq. 1) allocate(entropymat(nx,ny,neta))
        if(out_sel%bulk .eq. 1) allocate(bulkmat(nx,ny,neta))
        if(out_sel%pitt .eq. 1) allocate(ttmat(nx,ny,neta)) 
        if(out_sel%pitx .eq. 1) allocate(txmat(nx,ny,neta)) 
        if(out_sel%pity .eq. 1) allocate(tymat(nx,ny,neta)) 
        if(out_sel%pitz .eq. 1) allocate(tzmat(nx,ny,neta)) 
        if(out_sel%pixy .eq. 1) allocate(xymat(nx,ny,neta)) 
        if(out_sel%pixz .eq. 1) allocate(xzmat(nx,ny,neta)) 
        if(out_sel%piyz .eq. 1) allocate(yzmat(nx,ny,neta)) 
        if(out_sel%pixx .eq. 1) allocate(xxmat(nx,ny,neta)) 
        if(out_sel%piyy .eq. 1) allocate(yymat(nx,ny,neta)) 
        if(out_sel%pizz .eq. 1) allocate(zzmat(nx,ny,neta)) 
        if(out_sel%v0 .eq. 1) allocate(v0mat(nx,ny,neta)) 
        if(out_sel%dutdt .eq. 1) allocate(dutdtmat(nx,ny,neta))
        if(out_sel%duxdx .eq. 1) allocate(duxdxmat(nx,ny,neta))
        if(out_sel%duydy .eq. 1) allocate(duydymat(nx,ny,neta))
        if(out_sel%duzdz .eq. 1) allocate(duzdzmat(nx,ny,neta))
        if(out_sel%theta .eq. 1) allocate(thetamat(nx,ny,neta))
      else if(out_precision .eq. 4) then
        if(out_sel%density .eq. 1) allocate(rhomat4(nx,ny,neta))
        if(out_sel%vx .eq. 1) allocate(vxmat4(nx,ny,neta))
        if(out_sel%vy .eq. 1) allocate(vymat4(nx,ny,neta))
        if(out_sel%vz .eq. 1) allocate(vetamat4(nx,ny,neta))
        if(out_sel%pressure .eq. 1) allocate(prmat4(nx,ny,neta))
        if(out_sel%energy_density .eq. 1) allocate(epsmat4(nx,ny,neta))
        if(out_sel%temperature .eq. 1) allocate(tempmat4(nx,ny,neta))
        if(out_sel%entropy_density .eq. 1) allocate(entropymat4(nx,ny,neta))
        if(out_sel%bulk .eq. 1) allocate(bulkmat4(nx,ny,neta))
        if(out_sel%pitt .eq. 1) allocate(ttmat4(nx,ny,neta)) 
        if(out_sel%pitx .eq. 1) allocate(txmat4(nx,ny,neta)) 
        if(out_sel%pity .eq. 1) allocate(tymat4(nx,ny,neta)) 
        if(out_sel%pitz .eq. 1) allocate(tzmat4(nx,ny,neta)) 
        if(out_sel%pixy .eq. 1) allocate(xymat4(nx,ny,neta)) 
        if(out_sel%pixz .eq. 1) allocate(xzmat4(nx,ny,neta)) 
        if(out_sel%piyz .eq. 1) allocate(yzmat4(nx,ny,neta)) 
        if(out_sel%pixx .eq. 1) allocate(xxmat4(nx,ny,neta)) 
        if(out_sel%piyy .eq. 1) allocate(yymat4(nx,ny,neta)) 
        if(out_sel%pizz .eq. 1) allocate(zzmat4(nx,ny,neta)) 
        if(out_sel%v0 .eq. 1) allocate(v0mat4(nx,ny,neta)) 
        if(out_sel%dutdt .eq. 1) allocate(dutdtmat4(nx,ny,neta))
        if(out_sel%duxdx .eq. 1) allocate(duxdxmat4(nx,ny,neta))
        if(out_sel%duydy .eq. 1) allocate(duydymat4(nx,ny,neta))
        if(out_sel%duzdz .eq. 1) allocate(duzdzmat4(nx,ny,neta))
        if(out_sel%theta .eq. 1) allocate(thetamat4(nx,ny,neta))
      else
        write(*,*) "Sorry, I am not able to determine if output files have 4 or 8 bytes precision..."
        call exit(3)
      end if
 
      
      open (16, file=trim(adjustl(outdir))//'/tau.dat', status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
        write(*,*) trim(adjustl(outdir))//"/tau.dat cannot be created, so I'm forced to quit!"
        call exit(4)
      end if
      
      
      do it=ntstart,ntend
        nv3=it/1000+48
        nv2=it/100-it/1000*10+48
        nv1=it/10-it/100*10+48
        nv0=it-it/10*10+48
        no=CHAR(nv3)//CHAR(nv2)//CHAR(nv1)//CHAR(nv0)
        
        open (1, file='out'//no//'.dat', status='OLD', iostat=filerror, form='unformatted', access='stream')
        if (filerror .ne. 0) then
          write(*,*) 'out'//no//".dat cannot be opened, so I'm forced to quit!"
          close(1)
          call exit(5)
        end if

        read(1) time
        write(6,*)' time=',time

        write(16,*) it,time
        
      !read the variables
       if(out_precision .eq. 8) then
          do ieta=1,neta
            do iy=1,ny
              do ix=1,nx
                  if(out_sel%density .eq. 1) read(1) rhomat(ix,iy,ieta)
                  if(out_sel%vx .eq. 1) read(1) vxmat(ix,iy,ieta)
                  if(out_sel%vy .eq. 1) read(1) vymat(ix,iy,ieta)
                  if(out_sel%vz .eq. 1) read(1) vetamat(ix,iy,ieta)
                  if(out_sel%pressure .eq. 1) read(1) prmat(ix,iy,ieta)
                  if(out_sel%energy_density .eq. 1) read(1) epsmat(ix,iy,ieta)
                  if(out_sel%temperature .eq. 1) read(1) tempmat(ix,iy,ieta)
                  if(out_sel%entropy_density .eq. 1) read(1) entropymat(ix,iy,ieta)
                  if(out_sel%bulk .eq. 1) read(1) bulkmat(ix,iy,ieta)
                  if(out_sel%pitt .eq. 1) read(1) ttmat(ix,iy,ieta)
                  if(out_sel%pitx .eq. 1) read(1) txmat(ix,iy,ieta) 
                  if(out_sel%pity .eq. 1) read(1) tymat(ix,iy,ieta)
                  if(out_sel%pitz .eq. 1) read(1) tzmat(ix,iy,ieta) 
                  if(out_sel%pixy .eq. 1) read(1) xymat(ix,iy,ieta)
                  if(out_sel%pixz .eq. 1) read(1) xzmat(ix,iy,ieta) 
                  if(out_sel%piyz .eq. 1) read(1) yzmat(ix,iy,ieta)
                  if(out_sel%pixx .eq. 1) read(1) xxmat(ix,iy,ieta) 
                  if(out_sel%piyy .eq. 1) read(1) yymat(ix,iy,ieta)
                  if(out_sel%pizz .eq. 1) read(1) zzmat(ix,iy,ieta)
                  if(out_sel%v0 .eq. 1) read(1) v0mat(ix,iy,ieta) 
                  if(out_sel%dutdt .eq. 1) read(1) dutdtmat(ix,iy,ieta)
                  if(out_sel%duxdx .eq. 1) read(1) duxdxmat(ix,iy,ieta)
                  if(out_sel%duydy .eq. 1) read(1) duydymat(ix,iy,ieta)
                  if(out_sel%duzdz .eq. 1) read(1) duzdzmat(ix,iy,ieta)
                  if(out_sel%theta .eq. 1) read(1) thetamat(ix,iy,ieta)
              end do
            end do
          end do
       else
          do ieta=1,neta
            do iy=1,ny
              do ix=1,nx
                  if(out_sel%density .eq. 1) read(1) rhomat4(ix,iy,ieta)
                  if(out_sel%vx .eq. 1) read(1) vxmat4(ix,iy,ieta)
                  if(out_sel%vy .eq. 1) read(1) vymat4(ix,iy,ieta)
                  if(out_sel%vz .eq. 1) read(1) vetamat4(ix,iy,ieta)
                  if(out_sel%pressure .eq. 1) read(1) prmat4(ix,iy,ieta)
                  if(out_sel%energy_density .eq. 1) read(1) epsmat4(ix,iy,ieta)
                  if(out_sel%temperature .eq. 1) read(1) tempmat4(ix,iy,ieta)
                  if(out_sel%entropy_density .eq. 1) read(1) entropymat4(ix,iy,ieta)
                  if(out_sel%bulk .eq. 1) read(1) bulkmat4(ix,iy,ieta)
                  if(out_sel%pitt .eq. 1) read(1) ttmat4(ix,iy,ieta)
                  if(out_sel%pitx .eq. 1) read(1) txmat4(ix,iy,ieta) 
                  if(out_sel%pity .eq. 1) read(1) tymat4(ix,iy,ieta)
                  if(out_sel%pitz .eq. 1) read(1) tzmat4(ix,iy,ieta) 
                  if(out_sel%pixy .eq. 1) read(1) xymat4(ix,iy,ieta)
                  if(out_sel%pixz .eq. 1) read(1) xzmat4(ix,iy,ieta) 
                  if(out_sel%piyz .eq. 1) read(1) yzmat4(ix,iy,ieta)
                  if(out_sel%pixx .eq. 1) read(1) xxmat4(ix,iy,ieta) 
                  if(out_sel%piyy .eq. 1) read(1) yymat4(ix,iy,ieta)
                  if(out_sel%pizz .eq. 1) read(1) zzmat4(ix,iy,ieta)
                  if(out_sel%v0 .eq. 1) read(1) v0mat4(ix,iy,ieta) 
                  if(out_sel%dutdt .eq. 1) read(1) dutdtmat4(ix,iy,ieta)
                  if(out_sel%duxdx .eq. 1) read(1) duxdxmat4(ix,iy,ieta)
                  if(out_sel%duydy .eq. 1) read(1) duydymat4(ix,iy,ieta)
                  if(out_sel%duzdz .eq. 1) read(1) duzdzmat4(ix,iy,ieta)
                  if(out_sel%theta .eq. 1) read(1) thetamat4(ix,iy,ieta)
              end do
            end do
          end do
        end if 
       if(out_sel%density .eq. 1) then
          write(rhofile,'(a3,a4,a4)') 'RHO',no,'.dat'
          open (7, file=trim(adjustl(outdir))//'/'//rhofile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//rhofile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) rhomat(ix,iy,ieta)
            end do
           end do      
          end do
          close(7)
        end if

       if(out_sel%vx .eq. 1) then
          write(vxfile,'(a2,a4,a4)') 'VX',no,'.dat'
          open (7, file=trim(adjustl(outdir))//'/'//vxfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//vxfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) vxmat(ix,iy,ieta)
            end do
           end do      
          end do
          close(7)
        end if

       if(out_sel%vy .eq. 1) then
          write(vyfile,'(a2,a4,a4)') 'VY',no,'.dat'
          open (7, file=trim(adjustl(outdir))//'/'//vyfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//vyfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) vymat(ix,iy,ieta)
            end do
           end do      
          end do
          close(7)
        end if

       if(out_sel%vz .eq. 1) then
         write(vzfile,'(a2,a4,a4)') 'VZ',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//vzfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//vzfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) vetamat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pressure .eq. 1) then
         write(prfile,'(a2,a4,a4)') 'PR',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//prfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//prfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) prmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%energy_density .eq. 1) then
         write(epsfile,'(a3,a4,a4)') 'EPS',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//epsfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//epsfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) epsmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%temperature .eq. 1) then
         write(tempfile,'(a1,a4,a4)') 'T',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//tempfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//tempfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) tempmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%entropy_density .eq. 1) then
         write(entropyfile,'(a1,a4,a4)') 'S',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//entropyfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//entropyfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) entropymat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%bulk .eq. 1) then
          write(bulkfile,'(a4,a4,a4)') 'bulk',no,'.dat'
          open (7, file=trim(adjustl(outdir))//'/'//bulkfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//bulkfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) bulkmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pitt .eq. 1) then
         write(ttfile,'(a2,a4,a4)') 'tt',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//ttfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//ttfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) ttmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pitx .eq. 1) then
         write(txfile,'(a2,a4,a4)') 'tx',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//txfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//txfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) txmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pity .eq. 1) then
         write(tyfile,'(a2,a4,a4)') 'ty',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//tyfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//tyfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) tymat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pitz .eq. 1) then
         write(tzfile,'(a2,a4,a4)') 'tz',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//tzfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//tzfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) tzmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pixy .eq. 1) then
         write(xyfile,'(a2,a4,a4)') 'xy',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//xyfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//xyfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) xymat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pixz .eq. 1) then
         write(xzfile,'(a2,a4,a4)') 'xz',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//xzfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//xzfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) xzmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%piyz .eq. 1) then
         write(yzfile,'(a2,a4,a4)') 'yz',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//yzfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//yzfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) yzmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pixx .eq. 1) then
         write(xxfile,'(a2,a4,a4)') 'xx',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//xxfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//xxfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) xxmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%piyy .eq. 1) then
         write(yyfile,'(a2,a4,a4)') 'yy',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//yyfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//yyfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) yymat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%pizz .eq. 1) then
         write(zzfile,'(a2,a4,a4)') 'zz',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//zzfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//zzfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) zzmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%v0 .eq. 1) then
         write(v0file,'(a2,a4,a4)') 'u0',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//v0file, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//v0file, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) v0mat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%dutdt .eq. 1) then
         write(dutdtfile,'(a5,a4,a4)') 'dutdt',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//dutdtfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//dutdtfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) dutdtmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%duxdx .eq. 1) then
         write(duxdxfile,'(a5,a4,a4)') 'duxdx',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//duxdxfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//duxdxfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) duxdxmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%duydy .eq. 1) then
         write(duydyfile,'(a5,a4,a4)') 'duydy',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//duydyfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//duydyfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) duydymat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%duzdz .eq. 1) then
         write(duzdzfile,'(a5,a4,a4)') 'duzdz',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//duzdzfile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//duzdzfile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) duzdzmat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if

       if(out_sel%theta .eq. 1) then
         write(thetafile,'(a5,a4,a4)') 'theta',no,'.dat'
         open (7, file=trim(adjustl(outdir))//'/'//thetafile, status='REPLACE', iostat=filerror)
          if (filerror .ne. 0) then
            write(*,*) trim(adjustl(outdir))//'/'//thetafile, "cannot be created, so I'm forced to quit!"
            call exit(7)
          end if
          do ieta=1,neta
           do iy=1,ny
            do ix=1,nx
              write(7,*) thetamat(ix,iy,ieta)
            end do
           end do
          end do
          close(7)
        end if
      end do
      
      open (11, file=trim(adjustl(outdir))//'/x.dat', status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
        write(*,*) trim(adjustl(outdir))//'/x.dat', "cannot be created, so I'm forced to quit!"
        call exit(4)
      end if
      open (41, file=trim(adjustl(outdir))//'/x-index.dat', status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
        write(*,*) trim(adjustl(outdir))//'/x-index.dat', "cannot be created, so I'm forced to quit!"
        call exit(4)
      end if
 
      open (12, file=trim(adjustl(outdir))//'/y.dat', status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
        write(*,*) trim(adjustl(outdir))//'/y.dat', "cannot be created, so I'm forced to quit!"
        call exit(4)
      end if
      open (42, file=trim(adjustl(outdir))//'/y-index.dat', status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
        write(*,*) trim(adjustl(outdir))//'/y-index.dat', "cannot be created, so I'm forced to quit!"
        call exit(4)
      end if

      open (13, file=trim(adjustl(outdir))//'/eta.dat', status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
        write(*,*) trim(adjustl(outdir))//'/eta.dat', "cannot be created, so I'm forced to quit!"
        call exit(4)
      end if
      open (43, file=trim(adjustl(outdir))//'/eta-index.dat', status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
        write(*,*) trim(adjustl(outdir))//'/eta-index.dat', "cannot be created, so I'm forced to quit!"
        call exit(4)
      end if

      open (14, file=trim(adjustl(outdir))//'/tauxyeta.dat', status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
        write(*,*) trim(adjustl(outdir))//'/tauxyeta.dat', "cannot be created, so I'm forced to quit!"
        call exit(4)
      end if

      open (15, file='grid.dat', status='OLD', iostat=filerror, form='formatted')
      if (filerror .ne. 0) then
        write(*,*) "grid.dat cannot be opened, so I'm forced to quit!"
        call exit(2)
      end if
      
      !we read again this values written at the beginning of the file
      read(15,*) nx, ny, neta
      
      allocate(xvector(nx),yvector(ny),etavector(neta))

      !do ix=1,nx
      read(15,*) xvector!(ix)
      !end do
      !do iy=1,ny
      read(15,*) yvector!(iy)
      !end do
      !do ieta=1,neta
      read(15,*) etavector!(ieta)
      !end do
      
      do ix=1,nx
       write(11,*) xvector(ix)
       write(41,*) ix, xvector(ix)
      end do
      
      do iy=1,ny
       write(12,*) yvector(iy)
       write(42,*) iy, yvector(iy)
      end do
      
      do ieta=1,neta
       write(13,*) etavector(ieta)
       write(43,*) ieta, etavector(ieta)
      end do
      
      write(14,*)taustep
      write(14,*)xstep,nx
      write(14,*)ystep,ny
      write(14,*)etastep,neta

      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)
      close(41)
      close(42)
      close(43)

      call system('cp grid.dat '//trim(adjustl(outdir)))

      end program
