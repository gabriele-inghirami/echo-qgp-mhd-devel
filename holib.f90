! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
! *                                                                           *
! *  File: holib.f90                                                          *
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
! *                                                                           *
! *  Contributors:                                                            *
! *  Minor modifications: Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)*
! *                                                                           *
! *****************************************************************************

module holib
!-- REC, DER high-order procedures (only explicit routines)

  implicit none

contains

!-- REC PROCEDURES


! *****************************************************************************

subroutine holib_rectvd2(f,gl,gr,n)
!-- 2nd order TVD reconstruction
!-- Different limiters may be specified (mm2, mc2, ...)

  integer,intent( in)                   :: n
  real(8)   ,intent( in),dimension(-1:n+2) :: f
  real(8)   ,intent(out),dimension( 0:n  ) :: gl,gr

  integer :: i
  real(8)    :: dp

! here we use the mm2 limiter, comment this block and uncomment the next one if you wish to use mc2
  do i=1,n
    dp=.5*mm2(f(i)-f(i-1),f(i+1)-f(i))
    gl(i  )=f(i)+dp
    gr(i-1)=f(i)-dp
  end do
  gl(0)=f(0)+.5*mm2(f(0)-f(-1),f(1)-f(0))
  gr(n)=f(n+1)-.5*mm2(f(n+1)-f(n),f(n+2)-f(n+1))

! Uncomment theses line and comment the lines above if you wish to use mc2 instead of mm2
!  do i=1,n
!    dp=.5*mc2(f(i)-f(i-1),f(i+1)-f(i))
!    gl(i  )=f(i)+dp
!    gr(i-1)=f(i)-dp
!  end do
!  gl(0)=f(0)+.5*mc2(f(0)-f(-1),f(1)-f(0))
!  gr(n)=f(n+1)-.5*mc2(f(n+1)-f(n),f(n+2)-f(n+1))


end subroutine holib_rectvd2

! *****************************************************************************

subroutine holib_recceno3(f,gl,gr,n)
!-- 3rd order CENO reconstruction
!-- Different limiters may be specified (mm2, mc2, ...)

  real(8),parameter,dimension(3,3) :: &
  c=reshape((/3,-10,15,-1,6,3,3,6,-1/)/8.,shape=(/3,3/))

  integer,intent( in)                   :: n
  real(8)   ,intent( in),dimension(-2:n+3) :: f
  real(8)   ,intent(out),dimension( 0:n  ) :: gl,gr

  integer              :: i
  real(8)   ,dimension(3) :: d
  real(8)                 :: f1,f2,f3,f4,f5,p

  do i=0,n

    f1=f(i-2)
    f2=f(i-1)
    f3=f(i  )
    f4=f(i+1)
    f5=f(i+2)
! here we use the mm2 limiter, comment the following line and uncomment the next one if you wish to use mc2
    p=f3+.5*mm2(f3-f2,f4-f3)
!   p=f3+.5*mc2(f3-f2,f4-f3)
    d(1)=c(1,1)*f1+c(2,1)*f2+c(3,1)*f3-p
    d(2)=c(1,2)*f2+c(2,2)*f3+c(3,2)*f4-p
    d(3)=c(1,3)*f3+c(2,3)*f4+c(3,3)*f5-p
    gl(i)=p+ceno3(d)

    f1=f(i+3)
    f2=f(i+2)
    f3=f(i+1)
    f4=f(i  )
    f5=f(i-1)
! here we use the mm2 limiter, comment the following line and uncomment the next one if you wish to use mc2
    p=f3+.5*mm2(f3-f2,f4-f3)
!   p=f3+.5*mc2(f3-f2,f4-f3)
    d(1)=c(1,1)*f1+c(2,1)*f2+c(3,1)*f3-p
    d(2)=c(1,2)*f2+c(2,2)*f3+c(3,2)*f4-p
    d(3)=c(1,3)*f3+c(2,3)*f4+c(3,3)*f5-p
    gr(i)=p+ceno3(d)

  end do

end subroutine holib_recceno3

! *****************************************************************************

subroutine holib_recweno3(f,gl,gr,n)
!-- 3rd order WENO reconstruction

  real(8),parameter,dimension(2)   :: d=(/1,2/)/3.
  real(8),parameter,dimension(2,2) :: c=reshape((/-1,3,1,1/)/2.,shape=(/2,2/))
  real(8),parameter                :: eps=1.e-6

  integer,intent( in)                   :: n
  real(8)   ,intent( in),dimension(-1:n+2) :: f
  real(8)   ,intent(out),dimension( 0:n  ) :: gl,gr

  integer :: i
  real(8)    :: f1,f2,f3
  real(8)    :: alpha1,beta1,p1
  real(8)    :: alpha2,beta2,p2

!DEC$ VECTOR ALWAYS
  do i=0,n

    f1=f(i-1)
    f2=f(i  )
    f3=f(i+1)
    beta1=(f2-f1)**2
    beta2=(f3-f2)**2
    alpha1=d(1)/(beta1+eps)**2
    alpha2=d(2)/(beta2+eps)**2
    p1=c(1,1)*f1+c(2,1)*f2
    p2=c(1,2)*f2+c(2,2)*f3
    gl(i)=(alpha1*p1+alpha2*p2)/(alpha1+alpha2)

    f1=f(i+2)
    f2=f(i+1)
    f3=f(i  )
    beta1=(f2-f1)**2
    beta2=(f3-f2)**2
    alpha1=d(1)/(beta1+eps)**2
    alpha2=d(2)/(beta2+eps)**2
    p1=c(1,1)*f1+c(2,1)*f2
    p2=c(1,2)*f2+c(2,2)*f3
    gr(i)=(alpha1*p1+alpha2*p2)/(alpha1+alpha2)

  end do

end subroutine holib_recweno3

! *****************************************************************************

subroutine holib_recweno5(f,gl,gr,n)
!-- 5th order WENO reconstruction

  real(8),parameter,dimension(3)   :: d=(/1,10,5/)/16.
  real(8),parameter,dimension(3,3) :: &
    c=reshape((/3,-10,15,-1,6,3,3,6,-1/)/8.,shape=(/3,3/))
  real(8),parameter                :: eps=1.e-6

  integer,intent( in)                   :: n
  real(8)   ,intent( in),dimension(-2:n+3) :: f
  real(8)   ,intent(out),dimension( 0:n  ) :: gl,gr

  integer :: i
  real(8)    :: f1,f2,f3,f4,f5
  real(8)    :: alpha1,beta1,p1
  real(8)    :: alpha2,beta2,p2
  real(8)    :: alpha3,beta3,p3

!DEC$ VECTOR ALWAYS
  do i=0,n

    f1=f(i-2)
    f2=f(i-1)
    f3=f(i  )
    f4=f(i+1)
    f5=f(i+2)
    beta1=13.*(f3-2*f2+f1)**2+3.*(3.*f3-4.*f2+f1)**2
    beta2=13.*(f4-2*f3+f2)**2+3.*(f4-f2)**2
    beta3=13.*(f5-2*f4+f3)**2+3.*(f5-4.*f4+3.*f3)**2
    alpha1=d(1)/(beta1/12.+eps)**2
    alpha2=d(2)/(beta2/12.+eps)**2
    alpha3=d(3)/(beta3/12.+eps)**2
    p1=c(1,1)*f1+c(2,1)*f2+c(3,1)*f3
    p2=c(1,2)*f2+c(2,2)*f3+c(3,2)*f4
    p3=c(1,3)*f3+c(2,3)*f4+c(3,3)*f5
    gl(i)=(alpha1*p1+alpha2*p2+alpha3*p3)/(alpha1+alpha2+alpha3)

    f1=f(i+3)
    f2=f(i+2)
    f3=f(i+1)
    f4=f(i  )
    f5=f(i-1)
    beta1=13.*(f3-2*f2+f1)**2+3.*(3.*f3-4.*f2+f1)**2
    beta2=13.*(f4-2*f3+f2)**2+3.*(f4-f2)**2
    beta3=13.*(f5-2*f4+f3)**2+3.*(f5-4.*f4+3.*f3)**2
    alpha1=d(1)/(beta1/12.+eps)**2
    alpha2=d(2)/(beta2/12.+eps)**2
    alpha3=d(3)/(beta3/12.+eps)**2
    p1=c(1,1)*f1+c(2,1)*f2+c(3,1)*f3
    p2=c(1,2)*f2+c(2,2)*f3+c(3,2)*f4
    p3=c(1,3)*f3+c(2,3)*f4+c(3,3)*f5
    gr(i)=(alpha1*p1+alpha2*p2+alpha3*p3)/(alpha1+alpha2+alpha3)

  end do

end subroutine holib_recweno5

! *****************************************************************************

subroutine holib_recppm4(f,gl,gr,n)
!-- 4th order PPM reconstruction
!-- Different limiters may be specified (mm2, mc2, ...)

  integer,intent( in)                   :: n
  real(8)   ,intent( in),dimension(-1:n+2) :: f
  real(8)   ,intent(out),dimension( 0:n  ) :: gl,gr

  integer :: i
  real(8)    :: f1,f2,f3,f4

!DEC$ VECTOR ALWAYS
  do i=0,n

    f1=f(i-1)
    f2=f(i  )
    f3=f(i+1)
    f4=f(i+2)
! here we use the mm2 limiter, comment the following line and uncomment the next one if you wish to use mc2
    gl(i)=f2+.5*(f3-f2)-.125*(mm2(f4-f3,f3-f2)-mm2(f3-f2,f2-f1))
!   gl(i)=f2+.5*(f3-f2)-.125*(mc2(f4-f3,f3-f2)-mc2(f3-f2,f2-f1))

    f1=f(i+2)
    f2=f(i+1)
    f3=f(i  )
    f4=f(i-1)
! here we use the mm2 limiter, comment the following line and uncomment the next one if you wish to use mc2
    gr(i)=f2+.5*(f3-f2)-.125*(mm2(f4-f3,f3-f2)-mm2(f3-f2,f2-f1))
!   gr(i)=f2+.5*(f3-f2)-.125*(mc2(f4-f3,f3-f2)-mc2(f3-f2,f2-f1))

  end do

end subroutine holib_recppm4

! *****************************************************************************

subroutine holib_recmpe3(f,gl,gr,n)
!-- 3rd order explicit reconstruction

  real(8),parameter,dimension(3) :: d=(/-1,6,3/)/8.

  integer,intent( in)                   :: n
  real(8)   ,intent( in),dimension(-2:n+3) :: f
  real(8)   ,intent(out),dimension( 0:n  ) :: gl,gr

  integer :: i

  do i=0,n
    gl(i)=d(1)*f(i-1)+d(2)*f(i  )+d(3)*f(i+1)
    gr(i)=d(1)*f(i+2)+d(2)*f(i+1)+d(3)*f(i  )
  end do

!-- MP5 filter
  do i=0,n
    gl(i)=mp5(gl(i),f(i-2),f(i-1),f(  i),f(i+1),f(i+2))
    gr(i)=mp5(gr(i),f(i+3),f(i+2),f(i+1),f(i  ),f(i-1))
  end do

end subroutine holib_recmpe3

! *****************************************************************************

subroutine holib_recmpe5(f,gl,gr,n)
!-- 5th order explicit reconstruction

  real(8),parameter,dimension(5) :: d=(/3,-20,90,60,-5/)/128.

  integer,intent( in)                   :: n
  real(8)   ,intent( in),dimension(-2:n+3) :: f
  real(8)   ,intent(out),dimension( 0:n  ) :: gl,gr

  integer :: i

  do i=0,n
    gl(i)=d(1)*f(i-2)+d(2)*f(i-1)+d(3)*f(i  )+d(4)*f(i+1)+d(5)*f(i+2)
    gr(i)=d(1)*f(i+3)+d(2)*f(i+2)+d(3)*f(i+1)+d(4)*f(i  )+d(5)*f(i-1)
  end do

!-- MP5 filter
  do i=0,n
    gl(i)=mp5(gl(i),f(i-2),f(i-1),f(  i),f(i+1),f(i+2))
    gr(i)=mp5(gr(i),f(i+3),f(i+2),f(i+1),f(i  ),f(i-1))
  end do

end subroutine holib_recmpe5

! *****************************************************************************

subroutine holib_recmpe7(f,gl,gr,n)
!-- 7th order explicit reconstruction

  real(8),parameter,dimension(7) :: d=(/-5,42,-175,700,525,-70,7/)/1024.

  integer,intent( in)                   :: n
  real(8)   ,intent( in),dimension(-3:n+4) :: f
  real(8)   ,intent(out),dimension( 0:n  ) :: gl,gr

  integer :: i

  do i=0,n
    gl(i)=d(1)*f(i-3)+d(2)*f(i-2)+d(3)*f(i-1)+d(4)*f(i  ) &
         +d(5)*f(i+1)+d(6)*f(i+2)+d(7)*f(i+3)
    gr(i)=d(1)*f(i+4)+d(2)*f(i+3)+d(3)*f(i+2)+d(4)*f(i+1) &
         +d(5)*f(i  )+d(6)*f(i-1)+d(7)*f(i-2)
  end do

!-- MP5 filter
  do i=0,n
    gl(i)=mp5(gl(i),f(i-2),f(i-1),f(  i),f(i+1),f(i+2))
    gr(i)=mp5(gr(i),f(i+3),f(i+2),f(i+1),f(i  ),f(i-1))
  end do

end subroutine holib_recmpe7

! *****************************************************************************


!-- DER PROCEDURES


! *****************************************************************************

subroutine holib_dere2(f,g,n)
!-- 2nd order explicit derivation

  integer,intent( in)                :: n
  real(8),   intent( in),dimension(0:n) :: f
  real(8),   intent(out),dimension(1:n) :: g

  integer :: i

  do i=1,n
    g(i)=f(i)-f(i-1)
  end do

end subroutine holib_dere2

! *****************************************************************************

subroutine holib_dere4(f,g,n)
!-- 4th order explicit derivation

  real(8),parameter,dimension(4) :: d=(/1.0,-27.0,27.0,-1.0/)/24.

  integer,intent( in)                   :: n
  real(8),   intent( in),dimension(-1:n+1) :: f
  real(8),   intent(out),dimension( 1:n  ) :: g

  integer :: i

  do i=1,n
    g(i)=d(1)*f(i-2)+d(2)*f(i-1)+d(3)*f(i)+d(4)*f(i+1)
  end do

end subroutine holib_dere4

! *****************************************************************************

subroutine holib_dere6(f,g,n)
!-- 6th order explicit derivation

  real(8),parameter,dimension(6) :: d=(/-9,125,-2250,2250,-125,9/)/1920.

  integer,intent( in)                   :: n
  real(8),   intent( in),dimension(-2:n+2) :: f
  real(8),   intent(out),dimension( 1:n  ) :: g

  integer :: i

  do i=1,n
    g(i)=d(1)*f(i-3)+d(2)*f(i-2)+d(3)*f(i-1) &
        +d(4)*f(i  )+d(5)*f(i+1)+d(6)*f(i+2)
  end do

end subroutine holib_dere6

! *****************************************************************************

subroutine holib_dere8(f,g,n)
!-- 8th order explicit derivation

  real(8),parameter,dimension(8) :: d=(/75,-1029,8575,-128625, &
                                     128625,-8575,1029,-75/)/107520.

  integer,intent( in)                   :: n
  real(8),   intent( in),dimension(-3:n+3) :: f
  real(8),   intent(out),dimension( 1:n  ) :: g

  integer :: i

  do i=1,n
    g(i)=d(1)*f(i-4)+d(2)*f(i-3)+d(3)*f(i-2)+d(4)*f(i-1) &
        +d(5)*f(i  )+d(6)*f(i+1)+d(7)*f(i+2)+d(8)*f(i+3)
  end do

end subroutine holib_dere8

! *****************************************************************************

!-- LIMITERS AND FILTERS


! *****************************************************************************

function mm2(d1,d2)

  real(8)            :: mm2
  real(8),intent(in) :: d1,d2

  real(8) :: s1,s2

  s1=sign(real(1.,8),d1)
  s2=sign(real(1.,8),d2)

  mm2=.5*(s1+s2)*min(abs(d1),abs(d2))

end function mm2

! *****************************************************************************

function mc2(d1,d2)

  real(8)            :: mc2
  real(8),intent(in) :: d1,d2

  real(8) :: s1,s2

  s1=sign(real(1.,8),d1)
  s2=sign(real(1.,8),d2)

  mc2=.5*(s1+s2)*min(2*abs(d1),2*abs(d2),.5*abs(d1+d2))

end function mc2

! *****************************************************************************

function mm4(d1,d2,d3,d4)

  real(8)            :: mm4
  real(8),intent(in) :: d1,d2,d3,d4

  real(8) :: s1,s2,s3,s4

  s1=sign(real(1.,8),d1)
  s2=sign(real(1.,8),d2)
  s3=sign(real(1.,8),d3)
  s4=sign(real(1.,8),d4)

  mm4=.125*(s1+s2)*abs((s1+s3)*(s1+s4))*min(abs(d1),abs(d2),abs(d3),abs(d4))

end function mm4

! *****************************************************************************

function ceno3(d)

  real(8),parameter :: alpha=0.7

  real(8)                         :: ceno3
  real(8),intent(in),dimension(3) :: d

  integer,dimension(1) :: m
  integer,dimension(3) :: is
  real(8)   ,dimension(3) :: w

  ceno3=0.
  is(1:3)=sign(real(1.,8),d(1:3))
  if (abs(sum(is(1:3)))==3) then
    w(1:3)=abs(d(1:3))
    w(2)=alpha*w(2)
    m=minloc(w(1:3))
    ceno3=d(m(1))
  end if

end function ceno3

! *****************************************************************************

function mp5(f,f1,f2,f3,f4,f5)

  real(8)            :: mp5
  real(8),intent(in) :: f,f1,f2,f3,f4,f5

  real(8) :: ful,fmp,dm,d0,dp,fmd,flc,fmin,fmax

  mp5=f

  ful=f3+2*(f3-f2)
  fmp=f3+mm2(f4-f3,ful-f3)
  if ((f-f3)*(f-fmp)>0.) then
    dm=f1-2*f2+f3
    d0=f2-2*f3+f4
    dp=f3-2*f4+f5
    fmd=.5*(f3+f4)-.5*mm4(4*d0-dp,4*dp-d0,d0,dp)
    flc=f3+.5*(f3-f2)+(4./3.)*mm4(4*d0-dm,4*dm-d0,d0,dm)
    fmin=max(min(f3,f4,fmd),min(f3,ful,flc))
    fmax=min(max(f3,f4,fmd),max(f3,ful,flc))
    mp5=f+mm2(fmin-f,fmax-f)
  end if

end function mp5

! *****************************************************************************

end module holib

! *****************************************************************************
