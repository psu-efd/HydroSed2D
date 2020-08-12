c Aux functions 
c      distance between two points
       function dist2points(A, B)
          implicit none
          real*8 A(3), B(3), dist2points
          
          dist2points = dsqrt((A(1)-B(1))**2+
     *                        (A(2)-B(2))**2+
     *                        (A(3)-B(3))**2)
          return     
       end
