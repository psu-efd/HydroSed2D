!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    HydroSed2D Copyright (C) 2008 Xiaofeng Liu
!
!    License
!
!    This file is part of HydroSed2D.
!
!    HydroSed2D is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    HydroSed2D is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with HydroSed2D.  If not, see <http://www.gnu.org/licenses/>.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Aux subroutines 

!     vector dot product C = A & B
      subroutine vec_dot(A, B, C)
         implicit none
         real*8 A(3), B(3), C
         integer i

         C = 0.0
         do i = 1, 3
            C = C + A(i)*B(i)
         enddo
         return
      end

!     vector cross product C = AXB
      subroutine vec_cross(A, B, C)
         implicit none
         real*8 A(3), B(3), C(3)
 
         C(1) = A(2)*B(3) - A(3)*B(2)
         C(2) = A(3)*B(1) - A(1)*B(3)
         C(3) = A(1)*B(2) - A(2)*B(1)
         return
      end

!     vector normaliztion
      subroutine vec_normalize(A)
         implicit none
         real*8 A(3), temp

         call vec_mag(A, temp)
         if(dabs(temp).lt.1e-16) then 
            write(*,*) 'Try to normalize zero vector!'
            stop
         end if

         A(1) = A(1)/temp
         A(2) = A(2)/temp
         A(3) = A(3)/temp

         return
      end

!     vector magnitude
      subroutine vec_mag(A, C)
         implicit none
         real*8 A(3), C

         call vec_dot(A, A, C)
         C = dsqrt(C)
         
         return
       end
       

!      ascending order of 3 integer numbers
       subroutine orderThreeInt(A, B, C)
          implicit none
          integer*4 A, B, C, temp1, temp2, temp3

          temp1 = min(A, B, C)
          temp2 = max(A, B, C)
          temp3 = A + B + C - temp1 - temp2

          A = temp1
          B = temp3
          C = temp2          
       end


