c Aux functions 

c     vector dot product C = A & B
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

c     vector cross product C = AXB
      subroutine vec_cross(A, B, C)
         implicit none
         real*8 A(3), B(3), C(3)
 
         C(1) = A(2)*B(3) - A(3)*B(2)
         C(2) = A(3)*B(1) - A(1)*B(3)
         C(3) = A(1)*B(2) - A(2)*B(1)
         return
      end

c     vector normaliztion
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

c     vector magnitude
      subroutine vec_mag(A, C)
         implicit none
         real*8 A(3), C

         call vec_dot(A, A, C)
         C = dsqrt(C)
         
         return
       end
       
c      triangle face area
       subroutine triangle_area(A, B, C, area)
          implicit none
          real*8 A(3), B(3), C(3), area

          !not implemented yet
          area = 0.0
          return
       end

c      ascending order of N integer numbers
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


