! *****************************************************************************
! *                                                                           * 
! *  ECHO-QGP                                                                 * 
! *                                                                           * 
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015,2016 The ECHO-QGP team                                * 
! *                                                                           *
! *  File: glaubermc.f90                                                      *
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
! *  Authors: Andrea Beraudo (beraudo@to.infn.it)                             *
! *                                                                           *
! *  Contributors: Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)       *                                            
! *                                                                           *
! *  Acknowledgments:                                                         *
! *                                                                           *
! *****************************************************************************

module glaubermc
  use parallel
  use common
  implicit none
  integer, parameter :: ntrials=10000000
  real(8), parameter :: rmax=12.6d0 ! max radius of nucleons in fm
  real(8), parameter :: rmaxd=16.d0 ! max radius of nucleons in deuteron in fm
  real(8), parameter :: bmaxAA=20.d0 ! max b impact parameter in fm for AA collisions
  real(8), parameter :: bmaxdA=20.d0 ! max b impact parameter in fm for dA collisions
  real(8), parameter :: bmaxpA=13.d0 ! max b impact parameter in fm for pA collisions
  integer :: nevent,ncoll,npart,npartA,npartB
  real(8) ybeam

  
  real(8) xpartA(300),ypartA(300) !they will be not more than nA, but why to avoid a bit larger sizes?
  real(8) xspecA(300),yspecA(300)
  integer part_typeA(300), part_typeB(300), spec_typeA(300), spec_typeB(300)
  !for pA and dA 1 or 2 would be fine, but since they are not large arrays, let's avoid specific sizes for these cases
  real(8) xpartB(300),ypartB(300)
  real(8) xspecB(300),yspecB(300)
  real(8) xcoll(90000),ycoll(90000) !they will be not more than nA^2, but usually they will be at maximum a few thousands

  integer :: recl_partcoll_bin, recl_spectators_bin, recl_collisions_bin
  type event_data
    real(8) :: b
    integer :: ncoll
    integer :: npart
    integer :: npartA
    integer :: npartB
    integer :: nspecA
    integer :: nspecB
    integer :: index_partcoll_bin
    integer :: index_spectators_bin
    integer :: index_collisions_bin
  end type event_data

  type(event_data),dimension(:),allocatable :: events_data


contains


! **************************************************************************

subroutine generate_events

! subroutines to generate events for Glauber-MC initial conditions
! IMPORTANT: please, note that this subroutine usually is executed only by process 0 (so, if you use mpi, you have to broadcast)
      integer iseed
      integer previous_seed
      integer filerror, allocate_result,read_status
      !change it from -1 to the random seed of a previous run to reproduce its initial conditions

      integer i,j,iA,iB,ii,iimp,iconf !just counters DDD check if all of them are used

      real(8), allocatable, dimension(:,:) :: xA, yA, zA, xB, yB, zB
      real(8), allocatable, dimension(:,:) :: xcollAB, ycollAB, dsquare
      real(8), allocatable, dimension(:) :: xAcoll, xBcoll, dndnpart ! coordinates displaced of +/- b/2

      integer, allocatable, dimension(:) ::  ipartA ,ipartB ! flag to label participants
      integer, allocatable, dimension(:,:) :: icollAB ! flag to label collisions

      integer nA,nA1,nA2
      real(8) fmax
      real(8) b,r,rd,r1,r2,xcos,rcos,theta,phi,rphi,xcmA,xcmB,ycmA,ycmB,sumxA,sumyA,sumxB,sumyB,bmax,phin,pxB,pyB

      !for B field investigations
      integer nspecA, nspecB !nucleon spectators in nuclei A and B
      integer, allocatable, dimension(:) :: ispecA, ispecB !spectator nucleons
      integer :: record_length_real, record_length_integer
      integer :: idx_part, idx_spec, idx_coll, idp, ids, idc

      inquire (iolength=record_length_real) fmax
      inquire (iolength=record_length_integer) nA
      !
      recl_partcoll_bin=2*record_length_real+record_length_integer
      recl_spectators_bin=2*record_length_real+record_length_integer
      recl_collisions_bin=2*record_length_real
      
    
      open(unit=14,status='replace',file='summary.dat',iostat=filerror)
      if(filerror .ne. 0) then
          print*, 'Impossible to create file summary.dat, leaving...'
      end if

      open(unit=16,status='replace',file='partcoll.dat',iostat=filerror)      
      if(filerror .ne. 0) then
          print*, 'Impossible to create file partcoll.dat, leaving...'
      end if

      open(unit=17,status='replace',file='partcoll.bin',form='unformatted',access='direct',&
          &recl=recl_partcoll_bin,iostat=filerror)      
      if(filerror .ne. 0) then
          print*, 'Impossible to create file partcoll.bin, leaving...'
      end if

      open(unit=18,status='OLD',file='random_seed.dat',iostat=filerror)
      if (filerror .ne. 0) then
        call system_clock(iseed)
        call srand(iseed)
        open(unit=19,status='NEW',file='random_seed.dat',iostat=filerror)

        if(filerror .ne. 0) then
          print*, "ATTENTION, IT HAS NOT BEEN POSSIBLE TO STORE THE VALUE OF THE RANDOM SEED INTO THE FILE random_seed.dat"
        else
          print*,'New seed of the random number generator is:', iseed
          write(19,*) iseed
          close(19)
        end if
      else
        read(18,*) iseed
        close(18)
        print*,'Old seed of the random number generator was:', iseed
        call srand(iseed)
      end if

! added for B field study
      open(unit=20,status='replace',file='spectators.dat',iostat=filerror)
      if(filerror .ne. 0) then
          print*, 'Impossible to create file spectators.dat, leaving...'
      end if

      open(unit=21,status='replace',file='spectators.bin',form='unformatted',access='direct',&
          &recl=recl_spectators_bin,iostat=filerror)
      if(filerror .ne. 0) then
          print*, 'Impossible to create file spectators.bin, leaving...'
      end if

      open(unit=22,status='replace',file='collisions.bin',form='unformatted',access='direct',&
          &recl=recl_collisions_bin,iostat=filerror)
      if(filerror .ne. 0) then
          print*, 'Impossible to create file collisions.bin, leaving...'
      end if

! formats to read parameter file (to update if necessary...)
007	format(10x,f10.5)
008	format(10x,i1)
012	format(10x,i2)
014     format(17x,f8.4,9x,f7.5,13x,f8.6,28x,f8.6,46x,i3)
017     format(10x,a5)
018     format(10x,i5)
019     format(13x,a15)

009	format(a40,f10.3)
010	format(a40,3f12.3)
011	format(a40,e10.2)
112	format(a40,i10)
113     format(a40,a20)

      if(ienentr .ne. 0) then
         write(*,*) 'Careful!!! This initialization assumes that ienentr&
     &=0!! Please, stop the program, correct the parameters (usually, param.dat) file and&
     &relaunch the simulation again...'
      end if

      nA=int(projmass)
      nA1=nA
      select case (kind_of_collision)
        case(1)
          nA2=nA
          bmax=bmaxAA
        case(2)
          nA2=2
          bmax=bmaxdA
        case(3)
          nA2=1
          bmax=bmaxpA
        case default
          write(*,*) "Sorry, I was not able to understand which kind of collisions you wish to simulate..."
          write(*,*) "Please, check the parameter COLLISION into the parameters (usually, param.dat) file."
          call exit(1)
      end select

      if(kind_of_collision .eq. 3) then
        allocate(xA(1:nconf,1:nA1),yA(1:nconf,1:nA1),zA(1:nconf,1:nA1),dsquare(1:nA1,1:nA2),ipartA(1:nA1),dndnpart(1:nA1),&
               & stat=allocate_result)
        if(allocate_result /=0) then
          write(*,*)  "Proc.0 - Error, I can't allocate one or more of the arrays for Glauber-MonteCarlo initialization..."
          write(*,*)  "(source file glaubermc.f08)"
          call exit(1)
        end if
        xA=0.
        yA=0.
        zA=0.
        dndnpart=0
      else
        allocate(xA(1:nconf,1:nA1),yA(1:nconf,1:nA1),zA(1:nconf,1:nA1),xB(1:nconf,1:nA2),yB(1:nconf,1:nA2),zB(1:nconf,1:nA2),&
               & xcollAB(1:nA1,1:nA2), ycollAB(1:nA1,1:nA2),dsquare(1:nA1,1:nA2),ipartA(1:nA1),ipartB(1:nA2),icollAB(1:nA1,1:nA2),&
               & xAcoll(1:nA1),xBcoll(1:nA2), ispecA(1:nA1),ispecB(1:nA2),stat=allocate_result)
        if(allocate_result /=0) then
          write(*,*)  "Proc.0 - Error, I can't allocate one or more of the arrays for Glauber-MonteCarlo initialization..."
          write(*,*)  "(source file glaubermc.f08)"
          call exit(1)
        end if
        xA=0.
        yA=0.
        zA=0.
        xB=0.
        yB=0.
        zB=0.
      end if
      
 
      fmax=radius*radius*roze*2. ! to accept/reject
      nevent=0

      allocate(events_data(1:nconf),stat=allocate_result)
      if(allocate_result /=0) then
        write(*,*)  "Proc.0 - Error, I can't allocate the array of events_data structures..."
        write(*,*)  "(source file glaubermc.f08 for Glauber-MonteCarlo initialization)"
        call exit(1)
      end if


      ! we consider the first n=zelectrons as protons, the rest as neutrons
      do i=1,nconf !AAADDD the loop could include more lines of code instead of being repeated, left in this way for make easier
                   !direct comparison with previous versions of these subroutines
         iA=0
         ! extract radial position of nucleons from WS
         do j=1,ntrials
            r1=rand()
            r2=rand()
            r=r1*rmax  ! r1 within 0 and 1...
            if (woods(r,radius,delta,roze).le.(fmax*r2)) then
               continue
            else
               iA=iA+1
               ! evaluate x,y,z extracting angles
               rcos=rand()
               rphi=rand()
               xcos=2.d0*rcos-1.d0 
               theta=acos(xcos)
               phi=2.d0*pi*rphi

               xA(i,iA)=r*sin(theta)*cos(phi)
               yA(i,iA)=r*sin(theta)*sin(phi)
               zA(i,iA)=r*cos(theta)                        
       
            endif
            if (iA.eq.nA) exit !exits frm the j cycle and re-enters the i cycle with another nuclear configuration     
         enddo !close j=1,ntrials cycle
      end do  !close i=1,conf cycle

      select case(kind_of_collision)
        case(1) !AA collisions
          do i=1,nconf
          ! I extract then positions of nucleons in B
             iB=0
          ! extract radial position of nucleons from WS
             do j=1,ntrials
                r1=rand()
                r2=rand()
                r=r1*rmax  ! r1 within 0 and 1...
                if (woods(r,radius,delta,roze).le.(fmax*r2)) then
                   continue
                else
                   iB=iB+1
                   ! evaluate x,y,z extracting angles
                   rcos=rand()
                   rphi=rand()
                   xcos=2.d0*rcos-1.d0
                   theta=acos(xcos)
                   phi=2.d0*pi*rphi

                   xB(i,iB)=r*sin(theta)*cos(phi)
                   yB(i,iB)=r*sin(theta)*sin(phi)
                   zB(i,iB)=r*cos(theta)                        
  
                endif
                if (iB.eq.nA) exit      
             enddo !close j=1,ntrials cycle
          end do  !close i=1,conf cycle

        case(2) !dA collisions
          do i=1,nconf
          ! I extract then positions of nucleons in B
             iB=1 !we place the proton at the origin
             xB(i,iB)=0.
             yB(i,iB)=0.
             zB(i,iB)=0.
          ! extract radial position of nucleons from WS
             do j=1,ntrials
                r1=rand()
                r2=rand()
                rd=r1*rmaxd  ! r1 within 0 and 1...
                if (hulten(rd).le.(0.3*r2)) then
                   continue
                else
                   iB=iB+1
                   ! evaluate x,y,z extracting angles
                   rcos=rand()
                   rphi=rand()
                   xcos=2.d0*rcos-1.d0
                   theta=acos(xcos)
                   phi=2.d0*pi*rphi

                   xB(i,iB)=rd*sin(theta)*cos(phi)
                   yB(i,iB)=rd*sin(theta)*sin(phi)
                   zB(i,iB)=rd*cos(theta)                        
  
                endif
                if (iB.eq.nA2) exit      
             enddo !close j=1,ntrials cycle
          end do  !close i=1,conf cycle
      end select

      do i=1,nconf
         xcmA=0.
         ycmA=0.
         xcmB=0.
         ycmB=0.
         do j=1,nA1                  
            xcmA=xcmA+xA(i,j)
            ycmA=ycmA+yA(i,j)
         end do
         xcmA=xcmA/real(nA1,8)  ! center-of-mass coordinates
         ycmA=ycmA/real(nA1,8)
         if(kind_of_collision .ne. 3) then
           do j=1,nA2
              xcmB=xcmB+xB(i,j)
              ycmB=ycmB+yB(i,j)
           enddo
           xcmB=xcmB/real(nA2,8)
           ycmB=ycmB/real(nA2,8)         
         end if
         do j=1,nA1
            ! rewriting the coordinates in the cm frame                  
            xA(i,j)=xA(i,j)-xcmA
            yA(i,j)=yA(i,j)-ycmA
         end do
         if(kind_of_collision .ne. 3) then
           do j=1,nA2
              xB(i,j)=xB(i,j)-xcmB
              yB(i,j)=yB(i,j)-ycmB
              !write(12,*) xA(i,j),yA(i,j),xB(i,j),yB(i,j)
           enddo         
         end if
      end do  !close i=1,conf cycle
     
      nevent=0
      iconf=1 !iconf and nevent have the same value, they are just set at different times, but they are equivalent
      if(kind_of_collision .ne. 3) then
        do while(iconf .le. nconf)
          do iimp=1,nbcoll !start cycle nbcoll
           ! extract impact parameter
            r1=rand()
            r2=rand()
            b=r1*bmax  ! r1 within 0 and 1...
            if ((.not. fixed_b) .and. (b.le.(bmax*r2*1.001))) then ! linear function
              continue
            else
              if(fixed_b) b=bimpact
              ! shift the A and B positions along the x-axis
              do i=1,nA1
                xAcoll(i)=xA(iconf,i)+b/2.
              end do
              do i=1,nA2
                xBcoll(i)=xB(iconf,i)-b/2.                  
              enddo
            
              npartA=0  ! at each AA coll initialize npartA and npartB to 0
              npartB=0
              ncoll=0  ! at each AA coll initialize ncoll to 0
              ! at each AA, dA or pA collision initialize flags to 0
              ipartA=0 ! flag for participants
              ipartB=0
              icollAB=0 ! flag for collisions
          
              !protons spectators
              ispecA=0
              ispecB=0

              do iA=1,nA1
                 do iB=1,nA2
                    dsquare(iA,iB)=(xA(iconf,iA)-xB(iconf,iB)+b)**2+(yA(iconf,iA)-yB(iconf,iB))**2
                    if (dsquare(iA,iB).le.(sigma_in/pi)) then
                       ncoll=ncoll+1
                       ipartA(iA)=1
                       ipartB(iB)=1
                       icollAB(iA,iB)=1
                       xcollAB(iA,iB)=(xA(iconf,iA)+xB(iconf,iB))/2.
                       ycollAB(iA,iB)=(yA(iconf,iA)+yB(iconf,iB))/2.
                    endif
                 enddo              
              enddo

              nspecA=0
              do ii=1,nA1
                 npartA=npartA+ipartA(ii)
                 if(ipartA(ii) .eq. 0) then
                    nspecA=nspecA+1
                    ispecA(ii)=1
                  end if
              end do

              nspecB=0
              do ii=1,nA2
                 npartB=npartB+ipartB(ii)
                 if(ipartB(ii) .eq. 0) then
                    nspecB=nspecB+1
                    ispecB(ii)=1
                 end if
              enddo

              npart=npartA+npartB

              if ((ncoll.ge.1) .and. (npart .ge. min_participants)) then
                 nevent=nevent+1
                 write(14,*) 'ev=',nevent,' b=',b,' nc=',ncoll,' np=',npart, ' npA=',npartA, ' npB=',npartB, ' nspecA=',nspecA,&
                            &' nspecB=',nspecB

                 events_data(nevent)%b=b
                 events_data(nevent)%ncoll=ncoll
                 events_data(nevent)%npart=npart
                 events_data(nevent)%npartA=npartA
                 events_data(nevent)%npartB=npartB
                 events_data(nevent)%nspecA=nspecA
                 events_data(nevent)%nspecB=nspecB
                 if(nevent .eq. 1) then
                   events_data(nevent)%index_partcoll_bin=0
                   events_data(nevent)%index_spectators_bin=0
                   events_data(nevent)%index_collisions_bin=0
                 else
                   events_data(nevent)%index_partcoll_bin=events_data(nevent-1)%index_partcoll_bin+events_data(nevent)%npart
                   events_data(nevent)%index_spectators_bin=events_data(nevent-1)%index_spectators_bin+events_data(nevent)%nspecA&
                                                          &+events_data(nevent-1)%nspecB
                   events_data(nevent)%index_collisions_bin=events_data(nevent-1)%index_collisions_bin+events_data(nevent)%ncoll
                 end if

                 idx_part=events_data(nevent)%index_partcoll_bin
                 idx_spec=events_data(nevent)%index_spectators_bin
                 idx_coll=events_data(nevent)%index_collisions_bin


                 write(16,*) 'ev=',nevent,' b=',b,' nc=',ncoll,' np=',npart, ' npA=',npartA, ' npB=',npartB

                 if((nspecA .gt. 0) .or. (nspecB .gt. 0)) then
                   write(20,*) 'ev=',nevent,' b=', b, ' Aspectators=', nspecA,' Bspectators=', nspecB
                 end if

                 idp=0
                 ids=0
                 do iA=1,nA1
                   if (ipartA(iA).eq.1) then ! write participants of A
                     idp=idp+1
                     if(iA .le. zelectrons) then
                       write(16,*)'A',xAcoll(iA),yA(iconf,iA),'p'
                       write(17,rec=idx_part+idp) xAcoll(iA),yA(iconf,iA),0
                     else
                       write(16,*)'A',xAcoll(iA),yA(iconf,iA),'n'
                       write(17,rec=idx_part+idp) xAcoll(iA),yA(iconf,iA),1
                     end if
                   else if(ispecA(iA) .eq. 1) then
                       ids=ids+1
                       if(iA .le. zelectrons) then
                         write(20,*)'A',xAcoll(iA),yA(iconf,iA),'p'
                         write(21, rec=idx_spec+ids) xAcoll(iA),yA(iconf,iA),0
                       else
                         write(20,*)'A',xAcoll(iA),yA(iconf,iA),'n'
                         write(21, rec=idx_spec+ids) xAcoll(iA),yA(iconf,iA),1
                       end if
                   else
                     write(*,*) "Error found: a particle must be either a participant or a spectator..."
                     close(14)
                     close(16)
                     close(17)
                     close(20)
                     close(21)
                     close(22)
                     call exit(5)
                   end if
                 enddo

                 do iB=1,nA2
                    if (ipartB(iB).eq.1) then ! write participants of B
                     idp=idp+1
                     if(iB .le. zelectrons) then
                       write(16,*)'B',xBcoll(iB),yB(iconf,iB),'p'
                       write(17,rec=idx_part+idp) xBcoll(iB),yB(iconf,iB),0
                     else
                       write(16,*)'B',xBcoll(iB),yB(iconf,iB),'n'
                       write(17,rec=idx_part+idp) xBcoll(iB),yB(iconf,iB),1
                     end if
                    else if(ispecB(iB) .eq. 1) then
                      ids=ids+1
                      if(iB .le. zelectrons) then
                        write(20,*)'B',xBcoll(iB),yB(iconf,iB),'p'
                        write(21,rec=idx_spec+ids) xBcoll(iB),yB(iconf,iB),0
                      else
                        write(20,*)'B',xBcoll(iB),yB(iconf,iB),'n'
                        write(21,rec=idx_spec+ids) xBcoll(iB),yB(iconf,iB),1
                      end if
                    else
                     write(*,*) "Error found: a particle must be either a participant or a spectator..."
                     call exit(5)
                     close(14)
                     close(16)
                     close(17)
                     close(20)
                     close(21)
                     close(22)
                    endif
                 enddo

                 idc=0                 
                 do iA=1,nA1
                    do iB=1,nA2
                      if (icollAB(iA,iB).eq.1) then ! write collisions
                        idc=idc+1
                        write(16,*)'AB',xcollAB(iA,iB),ycollAB(iA,iB)
                        write(22,rec=idx_coll+idc) xcollAB(iA,iB),ycollAB(iA,iB)
                      endif
                    enddo
                 enddo

              else
                continue
              endif
            endif
          enddo ! end of the iimp=1,nbcoll cycle
          iconf=iconf+1
        enddo  ! end of the iconf=1,nconf cycle
      else  !begin case for pA collisions
        do iconf=1,nconf
          do iimp=1,nbcoll !start cycle nbcoll
           ! extract impact parameter
            r1=rand()
            r2=rand()
            b=r1*bmax  ! r1 within 0 and 1...
            if ((.not. fixed_b) .and. (b.le.(bmax*r2*1.02))) then ! linear function
              continue
            else
              if(fixed_b) b=bimpact
              phin=rand()
              phin=2.*pi*phin
              pxB=b*cos(phin)
              pyB=b*sin(phin)
            
              npartA=0  ! at each AA coll initialize npartA to 0
              ipartA=0 ! flag for participants
          
              do iA=1,nA1
                 do iB=1,nA2
                    dsquare(iA,iB)=(xA(iconf,iA)-pxB)**2+(yA(iconf,iA)-pyB)**2
                    if (dsquare(iA,iB).le.(sigma_in/pi)) then
                       ipartA(iA)=1
                    endif
                 enddo              
              enddo

              do ii=1,nA1
                 npartA=npartA+ipartA(ii)
              end do
              npart=npartA+1

              if (npart .ge. 2) then
                 nevent=nevent+1
                 write(14,*) 'ev=',nevent,' b=',b,' np=',npart
                 dndnpart(npart)=dndnpart(npart)+1

                 write(16,*) 'ev=',nevent,' b=',b,' np=',npart
                 idp=0
                 do iA=1,nA1
                   if (ipartA(iA).eq.1) then ! write participants of A
                    idp=idp+1
                    if(iA .lt. zelectrons) then
                      write(16,*)'A',xA(iconf,iA),yA(iconf,iA),"p"
                      write(17,rec=idx_part+idp) xA(iconf,iA),yA(iconf,iA),0
                    else
                      write(16,*)'A',xA(iconf,iA),yA(iconf,iA),"n"
                      write(17,rec=idx_part+idp) xA(iconf,iA),yA(iconf,iA),1
                    end if
                   endif
                 enddo

                 write(16,*)'p',pxB,pyB
                 write(22,rec=idx_coll+idc) pxB,pyB
              else
                continue
              endif
            endif
          enddo ! end of the iimp=1,nbcoll cycle
        enddo  ! end of the iconf=1,nconf cycle

      end if !end case for pA collisions
      

     close(14)
     close(16)
     close(20)
     close(17)
     close(21)
     close(22)

!    to uncomment if we want just to get the position of protons to compute an initial magnetic field
!     write(*,*) 'Events generated. Run finished.'
!     call exit(0)


end subroutine generate_events

! **************************************************************************

! Woods-Saxon distribution (no hard core!)
      real(8) function woods(r,radius,delta,roze)
      implicit none
      real(8) r,roze,radius,delta
      woods=r*r*roze/(exp((r-radius)/delta)+1.)
      return
      end

! **************************************************************************

! Hulten distribution (for deuteron)
      real(8) function hulten(r)
      implicit none
      real(8) r
      real(8), parameter :: a=0.228 !fm^-1
      real(8), parameter :: b=1.118 !fm^-1
      hulten=(exp(-a*r)-exp(-b*r))**2.
      return
      end

! **************************************************************************

      real(8) function feta(etas)
      implicit none
      real(8) etas
      real(8) half_deta

      half_deta=deta/2.

      if (abs(etas) .ge. ybeam) then
         feta=0.
         !the following is just an attemp to improve stability by avoiding jumps
         !however, please notice that without cutting the energy density at ybeam,
         !the profile designed in generate_energy has some issue
         !DDD feta=exp(-((abs(etas)-half_deta)**2.)/(2.*(sigeta**2.)))*(4.**(ybeam-abs(etas)))
      else
         if (abs(etas) .le. half_deta) then
            feta=1.
         else
            feta=exp(-((abs(etas)-half_deta)**2.)/(2.*(sigeta**2.)))
         end if
      end if
      return
      end
      
! **************************************************************************

subroutine generate_energy_density_profile(edens_array)

! program to plot the energy density
      use system_eqgp
      implicit none
     
      ! In the following:
      ! A: is the nucleus (usually Au), taken as right-moving
      ! B: is the second nucleus or the deuteron, taken as left-moving
      ! if one prefers the opposite kinematics, one has simply to invert the dependence on the rapidity

      !convention for type of particles (part_type and spec_type): 0=proton, 1=neutron 

      real(8), allocatable, dimension(:,:,:) :: edens_array
     
      character(3) event,impact,npart,AorB,numcoll,npartA,npartB

      real(8) b,epart,ecoll,edenscoll,edenspart,epartcoll,edenspartA,edenspartB

      real(8) xprot,yprot, ZZ_pA, AA_pA, ynucl, ycm, ecm, enucleus, eproton !variables used for the pA collisions case

      integer i,ievrecord,numpart,ix,iy,iz,iev,ncoll,numpartA,numpartB


      integer error_flag, allocate_result, opening_error

      real(8) cmx, cmy, cmx_p, cmx_c, cmy_p, cmy_c

      integer :: nucltype_bitbucket, index_partcoll_bin, index_collisions_bin, index_spectators_bin

      iev=id_of_the_run     
      if(pe0) then !oly the first processor knows the data of the events
        b=events_data(iev)%b
        bimpact=b  !it can change if the impact parameter is not fixed. We report and print this value in config_summary.dat
        ncoll=events_data(iev)%ncoll
        numpart=events_data(iev)%npart
        numpartA=events_data(iev)%npartA
        numpartB=events_data(iev)%npartB
        index_partcoll_bin=events_data(iev)%index_partcoll_bin
        index_collisions_bin=events_data(iev)%index_collisions_bin
      end if

      if(prl) then
        call MPI_Bcast(b,1,mpi_realtype,0,icomm,ierr)
        call MPI_Bcast(bimpact,1,mpi_realtype,0,icomm,ierr)
        call MPI_Bcast(ncoll,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(numpart,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(numpartA,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(numpartB,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(index_partcoll_bin,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(index_collisions_bin,1,mpi_integer,0,icomm,ierr)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end if 
 
      error_flag=0
      ybeam=log(rads/am)
      cmx=0.
      cmy=0.
      cmx_p=0.
      cmy_p=0.
      cmx_c=0.
      cmy_c=0.

      if(kind_of_collision .eq. 3) then
        yprot=log(2*eprot/am)
        AA_pA=real(projmass,8)
        ZZ_pA=real(zelectrons,8)
        ynucl=yprot+log(ZZ_pA/AA_pA)
        ycm=log(AA_pA/ZZ_pA)/2.d0
        ecm=2*eprot*sqrt(ZZ_pA/AA_pA)
      end if
        
    
  
      if(pe0) then !only the first processor manages the files
         open(unit=2,file='partcoll.bin',status='old',iostat=opening_error,access='direct',recl=recl_partcoll_bin)
         if(opening_error .ne. 0) then
            write(*,*) "Sorry, it was not possible to open the file partcoll.bin. Leaving..."
            call exit(1)
         end if

         open(unit=3,file='collisions.bin',status='old',iostat=opening_error,access='direct',recl=recl_collisions_bin)
         if(opening_error .ne. 0) then
            write(*,*) "Sorry, it was not possible to open the file partcoll.bin. Leaving..."
            call exit(1)
         end if

         write(*,*) "Reading event:", id_of_the_run

         if(kind_of_collision .le. 2) then
           do ievrecord=1,numpartA
              read(2,rec=index_partcoll_bin+ievrecord) xpartA(ievrecord),ypartA(ievrecord), nucltype_bitbucket
           enddo
           do ievrecord=1,numpartB
              read(2,rec=index_partcoll_bin+numpartA+ievrecord) xpartB(ievrecord),ypartB(ievrecord), nucltype_bitbucket
           enddo
           do ievrecord=1,ncoll
              read(3, rec=index_collisions_bin+ievrecord) xcoll(ievrecord),ycoll(ievrecord)
           end do
         else
           !read data for the pA case
           do ievrecord=1,numpart
              read(2,rec=index_partcoll_bin+ievrecord) xpartA(ievrecord),ypartA(ievrecord),nucltype_bitbucket
           enddo
         end if
      end if !end if pe0

  !processor 0 informs the other processors
  if(prl) then
    call MPI_Bcast(xpartA, numpartA, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(ypartA, numpartA, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(xpartB, numpartB, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(ypartB, numpartB, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(xcoll, ncoll, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(ycoll, ncoll, mpi_realtype, 0, icomm, ierr)
      
    !call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end if

      if(kind_of_collision .ne. 3) then
        do i=1,numpartA
           cmx_p=cmx_p+xpartA(i)
           cmy_p=cmy_p+ypartA(i)
        end do
        do i=1,numpartB
           cmx_p=cmx_p+xpartB(i)
           cmy_p=cmy_p+ypartB(i)
        end do
        do i=1,ncoll
           cmx_c=cmx_c+xcoll(i)
           cmy_c=cmy_c+ycoll(i)   
        end do
        cmx=((1-ah)*cmx_p+ah*cmx_c)/((1-ah)*(numpartA+numpartB)+ah*ncoll)
        cmy=((1-ah)*cmy_p+ah*cmy_c)/((1-ah)*(numpartA+numpartB)+ah*ncoll)

         do iz=iz1,iz2
            do iy=iy1,iy2
               do ix=ix1,ix2
                  edenspartA=0.
                  edenspartB=0.
                  do i=1,numpartA
                     edenspartA=edenspartA+exp(-((x(ix)-xpartA(i)+cmx)**2+(y(iy)-ypartA(i)+cmy)**2)/2./sig_mc**2)
                  enddo
                  do i=1,numpartB
                     edenspartB=edenspartB+exp(-((x(ix)-xpartB(i)+cmx)**2+(y(iy)-ypartB(i)+cmy)**2)/2./sig_mc**2)
                  enddo
                  edenspart=edenspartA*(ybeam+z(iz))/ybeam+edenspartB*(ybeam-z(iz))/ybeam
                  edenscoll=0.
                  do i=1,ncoll
                     edenscoll=edenscoll+exp(-((x(ix)-xcoll(i)+cmx)**2+(y(iy)-ycoll(i)+cmy)**2)/2./sig_mc**2)
                  enddo

                  edens_array(ix,iy,iz)=kappa/(2.d0*pi*sig_mc**2)*(ah*edenscoll+(1-ah)*edenspart)*feta(z(iz))
               enddo
            enddo
         end do
      else !case of pA collisions
         do i=1,numpartA
            cmx_p=cmx_p+xpartA(i)+xprot
            cmy_p=cmy_p+ypartA(i)+yprot
         end do
         cmx=cmx_p/(numpartA+1)
         cmy=cmy_p/(numpartA+1)
         do iz=iz1,iz2
            do iy=iy1,iy2
               do ix=ix1,ix2
                  enucleus=0.
                  do i=1,numpart-1
                     enucleus=enucleus+exp(-((x(ix)-xpartA(i)+cmx)**2+(y(iy)-ypartA(i)+cmy)**2)/2./sig_mc**2)
                  enddo
                  eproton=exp(-((x(ix)-xprot+cmx)**2+(y(iy)-yprot+cmy)**2)/2./sig_mc**2)
                  edenspart=eproton*(ybeam+z(iz))/ybeam+enucleus*(ybeam-z(iz))/ybeam
                  edens_array(ix,iy,iz)=edenspart*kappa/2./pi/sig_mc**2*feta(z(iz))

               enddo
            enddo
         end do
      !end of case of pA collisions
      end if 
      

      if(pe0) then
        close(2)
        close(3)
      end if

end subroutine generate_energy_density_profile

! **************************************************************************

subroutine generate_B_field_profile()

! program to plot the energy density
      use system_eqgp
      implicit none
     

      ! In the following:
      ! A: is the nucleus (usually Au), taken as right-moving
      ! B: is the second nucleus or the deuteron, taken as left-moving
      ! if one prefers the opposite kinematics, one has simply to invert the dependence on the rapidity

      !convention for type of particles (part_type and spec_type): 0=proton, 1=neutron 

      real(8), allocatable, dimension(:,:,:) :: edens_array
     
      character(3) event,impact,npart,AorB,numcoll,npartA,npartB,nspecA,nspecB

      real(8) b, tt, xcoord,ycoord,zcoord,zpos,vel, xt2, xt, D, sqrD, A, Bphi, Br, Bz, toMilne, toMilne_z, phi

      real(8) xprot,yprot, ZZ_pA, AA_pA, ynucl, ycm, ecm, enucleus, eproton !variables used for the pA collisions case

      integer i,ievrecord,numpart,ix,iy,iz,iev,ncoll,numpartA,numpartB,numspecA,numspecB


      integer error_flag, allocate_result, opening_error

      integer :: index_partcoll_bin, index_collisions_bin, index_spectators_bin

      iev=id_of_the_run     
      if(pe0) then !oly the first processor knows the data of the events
        b=events_data(iev)%b
        ncoll=events_data(iev)%ncoll
        numpart=events_data(iev)%npart
        numpartA=events_data(iev)%npartA
        numpartB=events_data(iev)%npartB
        numspecA=events_data(iev)%nspecA
        numspecB=events_data(iev)%nspecB
        index_partcoll_bin=events_data(iev)%index_partcoll_bin
        index_spectators_bin=events_data(iev)%index_spectators_bin
        index_collisions_bin=events_data(iev)%index_collisions_bin
      end if

      if(prl) then
        call MPI_Bcast(b,1,mpi_realtype,0,icomm,ierr)
        call MPI_Bcast(ncoll,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(numpart,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(numpartA,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(numpartB,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(numspecA,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(numspecB,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(index_partcoll_bin,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(index_spectators_bin,1,mpi_integer,0,icomm,ierr)
        call MPI_Bcast(index_collisions_bin,1,mpi_integer,0,icomm,ierr)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end if 

      error_flag=0
      ybeam=log(rads/am)
  
      if(pe0) then !only the first processor manages the files
         open(unit=2,file='partcoll.bin',status='old',iostat=opening_error,access='direct',recl=recl_partcoll_bin)
         if(opening_error .ne. 0) then
            write(*,*) "Sorry, it was not possible to open the file partcoll.bin. Leaving..."
            call exit(1)
         end if

         open(unit=3,file='collisions.bin',status='old',iostat=opening_error,access='direct',recl=recl_collisions_bin)
         if(opening_error .ne. 0) then
            write(*,*) "Sorry, it was not possible to open the file partcoll.bin. Leaving..."
            call exit(1)
         end if

         open(unit=4,file='spectators.bin',status='old',iostat=opening_error,access='direct',recl=recl_spectators_bin)
         if(opening_error .ne. 0) then
            write(*,*) "Sorry, it was not possible to open the file spectators.bin. Leaving..."
            call exit(1)
         end if

         write(*,*) "Reading event:", id_of_the_run

         if(kind_of_collision .le. 2) then
           do ievrecord=1,numpartA
              read(2,rec=index_partcoll_bin+ievrecord) xpartA(ievrecord),ypartA(ievrecord), part_typeA(ievrecord)
           enddo
           do ievrecord=1,numpartB
              read(2,rec=index_partcoll_bin+numpartA+ievrecord) xpartB(ievrecord),ypartB(ievrecord), part_typeB(ievrecord)
           enddo
           do ievrecord=1,ncoll
              read(3,rec= index_collisions_bin+ievrecord) xcoll(ievrecord),ycoll(ievrecord)
           end do
           do ievrecord=1,numspecA
              read(4,rec=index_spectators_bin+ievrecord) xspecA(ievrecord),yspecA(ievrecord), spec_typeA(ievrecord)
           enddo
           do ievrecord=1,numspecB
              read(4,rec=index_spectators_bin+numspecA+ievrecord) xspecB(ievrecord),yspecB(ievrecord), spec_typeB(ievrecord)
           enddo
         else
           !actually, this message should never be printed, because there is another check in init.f90 just before calling
           !this subroutine
           write(*,*) "Sorry, but the creation of B fields at runtime is not possible with pA collisions. I quit."
           call exit(1)
         end if

   end if !end pe0
  !processor 0 informs the other processors

  if(prl) then
    call MPI_Bcast(numpartA,1,mpi_integer,0,icomm,ierr)
    call MPI_Bcast(numpartB,1,mpi_integer,0,icomm,ierr)
    call MPI_Bcast(numspecA,1,mpi_integer,0,icomm,ierr)
    call MPI_Bcast(numspecB,1,mpi_integer,0,icomm,ierr)
    call MPI_Bcast(ncoll,1,mpi_integer,0,icomm,ierr)
    call MPI_Bcast(xpartA, numpartA, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(ypartA, numpartA, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(xspecA, numspecA, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(yspecA, numspecA, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(xpartB, numpartB, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(ypartB, numpartB, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(xspecB, numspecB, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(yspecB, numspecB, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(xcoll, ncoll, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(ycoll, ncoll, mpi_realtype, 0, icomm, ierr)
    call MPI_Bcast(spec_typeA, numspecA, mpi_integer, 0, icomm, ierr)
    call MPI_Bcast(spec_typeB, numspecB, mpi_integer, 0, icomm, ierr)
    call MPI_Bcast(part_typeA, numpartA, mpi_integer, 0, icomm, ierr)
    call MPI_Bcast(part_typeB, numpartB, mpi_integer, 0, icomm, ierr)
      
    !call MPI_Barrier(MPI_COMM_WORLD,ierr)

  end if

      vel=sqrt(1. - 1./(gamma_col**2))

      do iz=iz1,iz2
        if(coordinates .eq. 1) then
          zcoord=z(iz)
          tt=tmin
          toMilne=1
          toMilne_z=1
        else
          tt=tmin*cosh(z(iz))        
          zcoord=tmin*sinh(z(iz))        
          toMilne=cosh(z(iz))
          toMilne_z=tmin
        end if
        do iy=iy1,iy2
          do ix=ix1,ix2
             zpos=tt*vel+zcoord
             do ievrecord=1,numspecA
                if(spec_typeA(ievrecord) .eq. 0) then
                  xcoord=x(ix)-xspecA(ievrecord)
                  ycoord=y(iy)-yspecA(ievrecord)
                  phi=atan2(ycoord,xcoord)
                  xt2=xcoord**2+ycoord**2
                  xt=sqrt(xt2);
                  D=(gamma_col*zpos)**2+xt2;
                  sqrD=sqrt(D);
                  A=(sigma_el_cond*vel*gamma_col/2)*(gamma_col*zpos-sqrD)/hbar;

                  if(magfield_type .eq. 3) then
                   Bphi=0.d0
                  else
                   Bphi=B_amplify_factor*Qbig/(4*pi)*sqrh*(vel*gamma_col*xt/(D**1.5))*(1+sigma_el_cond*vel*gamma_col*sqrD/&
                       &(2*hbar))*exp(A)
                  end if
                  if(magfield_type .eq. 2) then
                   Br=0.d0
                   Bz=0.d0
                  else               
                   Br=B_amplify_factor*(-sigma_chiral_cond)*Qbig/(8*pi)*(vel*gamma_col**2*xt/(D**1.5))*(gamma_col*zpos+A*sqrD)&
                     &*exp(A)/sqrh
                   Bz=B_amplify_factor*(sigma_chiral_cond)*Qbig/(8*pi)*(vel*gamma_col/(D**1.5))*((gamma_col*zpos)**2*&
                     &(1+sigma_el_cond*vel*gamma_col*sqrD/(2*hbar))+D*(1-sigma_el_cond*vel*gamma_col*sqrD/(2*hbar)))*exp(A)/sqrh
                  end if
!                  Bphi=0.!GUUU
                  v(ix,iy,iz,kbx)=v(ix,iy,iz,kbx)+(Bphi*(-sin(phi))+Br*cos(phi))/toMilne
                  v(ix,iy,iz,kby)=v(ix,iy,iz,kby)+(Bphi*cos(phi)+Br*sin(phi))/toMilne
                  v(ix,iy,iz,kbz)=v(ix,iy,iz,kbz)+Bz/toMilne_z
                else !the protons are registered before neutrons, so if we find a neutron we're done and we can exit the cycle
                  exit
                end if
             end do

             do ievrecord=1,numpartA
                if(part_typeA(ievrecord) .eq. 0) then
                  xcoord=x(ix)-xpartA(ievrecord)
                  ycoord=y(iy)-ypartA(ievrecord)
                  phi=atan2(ycoord,xcoord)
                  xt2=xcoord**2+ycoord**2
                  xt=sqrt(xt2);
                  D=(gamma_col*zpos)**2+xt2;
                  sqrD=sqrt(D);
                  A=(sigma_el_cond*vel*gamma_col/2)*(gamma_col*zpos-sqrD)/hbar;

                  if(magfield_type .eq. 3) then
                   Bphi=0.d0
                  else
                   Bphi=B_amplify_factor*pw*Qbig/(4*pi)*sqrh*(vel*gamma_col*xt/(D**1.5))*(1+sigma_el_cond*vel*gamma_col*sqrD/&
                       &(2*hbar))*exp(A)
                  end if
                  if(magfield_type .eq. 2) then
                   Br=0.d0
                   Bz=0.d0
                  else               
                   Br=B_amplify_factor*pw*(-sigma_chiral_cond)*Qbig/(8*pi)*(vel*gamma_col**2*xt/(D**1.5))*(gamma_col*zpos+A*sqrD)&
                     &*exp(A)/sqrh
                   Bz=B_amplify_factor*pw*(sigma_chiral_cond)*Qbig/(8*pi)*(vel*gamma_col/(D**1.5))*((gamma_col*zpos)**2*&
                     &(1+sigma_el_cond*vel*gamma_col*sqrD/(2*hbar))+D*(1-sigma_el_cond*vel*gamma_col*sqrD/(2*hbar)))*exp(A)/sqrh
                  end if
!                  Bphi=0.!GUUU
                  v(ix,iy,iz,kbx)=v(ix,iy,iz,kbx)+(Bphi*(-sin(phi))+Br*cos(phi))/toMilne
                  v(ix,iy,iz,kby)=v(ix,iy,iz,kby)+(Bphi*cos(phi)+Br*sin(phi))/toMilne
                  v(ix,iy,iz,kbz)=v(ix,iy,iz,kbz)+Bz/toMilne_z
                else !the protons are registered before neutrons, so if we find a neutron we're done and we can exit the cycle
                  exit
                end if
             end do
             
             zpos=tt*vel-zcoord
             do ievrecord=1,numspecB
                if(spec_typeB(ievrecord) .eq. 0) then
                  xcoord=x(ix)-xspecB(ievrecord)
                  ycoord=y(iy)-yspecB(ievrecord)
                  phi=atan2(ycoord,xcoord)
                  xt2=xcoord**2+ycoord**2
                  xt=sqrt(xt2);
                  D=(gamma_col*zpos)**2+xt2;
                  sqrD=sqrt(D);
                  A=(sigma_el_cond*vel*gamma_col/2)*(gamma_col*zpos-sqrD)/hbar;

                  if(magfield_type .eq. 3) then
                   Bphi=0.d0
                  else
                   Bphi=B_amplify_factor*Qbig/(4*pi)*sqrh*(vel*gamma_col*xt/(D**1.5))*(1+sigma_el_cond*vel*gamma_col*sqrD/&
                       &(2*hbar))*exp(A)
                  end if
                  if(magfield_type .eq. 2) then
                   Br=0.d0
                   Bz=0.d0
                  else               
                   Br=B_amplify_factor*(-sigma_chiral_cond)*Qbig/(8*pi)*(vel*gamma_col**2*xt/(D**1.5))*(gamma_col*zpos+A*sqrD)&
                     &*exp(A)/sqrh
                   Bz=-B_amplify_factor*(sigma_chiral_cond)*Qbig/(8*pi)*(vel*gamma_col/(D**1.5))*((gamma_col*zpos)**2*&
                     &(1+sigma_el_cond*vel*gamma_col*sqrD/(2*hbar))+D*(1-sigma_el_cond*vel*gamma_col*sqrD/(2*hbar)))*exp(A)/sqrh
                  end if
!                  Bphi=0.!GUUU
                  v(ix,iy,iz,kbx)=v(ix,iy,iz,kbx)+(Bphi*sin(phi)+Br*cos(phi))/toMilne
                  v(ix,iy,iz,kby)=v(ix,iy,iz,kby)+(Bphi*(-cos(phi))+Br*sin(phi))/toMilne
                  v(ix,iy,iz,kbz)=v(ix,iy,iz,kbz)+Bz/toMilne_z
                else !the protons are registered before neutrons, so if we find a neutron we're done and we can exit the cycle
                  exit
                end if
             end do

             do ievrecord=1,numpartB
                if(part_typeB(ievrecord) .eq. 0) then
                  xcoord=x(ix)-xpartB(ievrecord)
                  ycoord=y(iy)-ypartB(ievrecord)
                  phi=atan2(ycoord,xcoord)
                  xt2=xcoord**2+ycoord**2
                  xt=sqrt(xt2);
                  D=(gamma_col*zpos)**2+xt2;
                  sqrD=sqrt(D);
                  A=(sigma_el_cond*vel*gamma_col/2)*(gamma_col*zpos-sqrD)/hbar;

                  if(magfield_type .eq. 3) then
                   Bphi=0.d0
                  else
                   Bphi=B_amplify_factor*pw*Qbig/(4*pi)*sqrh*(vel*gamma_col*xt/(D**1.5))*(1+sigma_el_cond*vel*gamma_col*sqrD/&
                       &(2*hbar))*exp(A)
                  end if
                  if(magfield_type .eq. 2) then
                   Br=0.d0
                   Bz=0.d0
                  else               
                   Br=B_amplify_factor*pw*(-sigma_chiral_cond)*Qbig/(8*pi)*(vel*gamma_col**2*xt/(D**1.5))*(gamma_col*zpos+A*sqrD)&
                     &*exp(A)/sqrh
                   Bz=-B_amplify_factor*pw*(sigma_chiral_cond)*Qbig/(8*pi)*(vel*gamma_col/(D**1.5))*((gamma_col*zpos)**2*&
                     &(1+sigma_el_cond*vel*gamma_col*sqrD/(2*hbar))+D*(1-sigma_el_cond*vel*gamma_col*sqrD/(2*hbar)))*exp(A)/sqrh
                  end if
!                  Bphi=0.!GUUU
                  v(ix,iy,iz,kbx)=v(ix,iy,iz,kbx)+(Bphi*sin(phi)+Br*cos(phi))/toMilne
                  v(ix,iy,iz,kby)=v(ix,iy,iz,kby)+(Bphi*(-cos(phi))+Br*sin(phi))/toMilne
                  v(ix,iy,iz,kbz)=v(ix,iy,iz,kbz)+Bz/toMilne_z
                else !the protons are registered before neutrons, so if we find a neutron we're done and we can exit the cycle
                  exit
                end if
             end do
          end do
        end do
      end do

      close(2)
      close(3)
      close(4)

end subroutine generate_B_field_profile

! **************************************************************************

end module
