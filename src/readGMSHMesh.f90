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

!   Read mesh from GMSH format files
	subroutine readGMSHMesh()
	USE COMMON_MODULE
	implicit none

	integer*4 i
	integer*4 nElementType
	integer*4 iTemp1, iTemp2
	character*100 string
    integer*4 iDummy, nElements
	integer*4 points_temp(3), mark_temp


     !read data files from GMSH format file of face from GMSH
	 if(platform.eq.1) then  !windows platform
		tempfilename=trim(trim(dirname) // '\input\' // meshfilename)
	 else                    !linux platform 
	    tempfilename=trim(trim(dirname) // '/input/' // meshfilename)
	 endif
       
	   open(1, file = tempfilename, status = 'old', iostat = ios)
       if ( ios /= 0 ) then
           ierror = 1
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'Could not open the input file'
           
		   stop
       end if

     !read out 4 line of comments
     read(1,*) string
	 read(1,*) string
	 read(1,*) string
	 read(1,*) string
       
	 read(1,*) nNodes

       do i = 1, nNodes
           read(1, *) iDummy,pcoor(i,1),pcoor(i,2),pcoor(i,3)
       enddo
	
	!read out 2 line of comments 
       read(1,*) string
       read(1,*) string

	!read elements: lines and triangles
	 nElementType=1
	 iTemp1=0
	 
	 read(1,*) nElements

	 do i=1, nElements
 		 read(1,*) iDummy,nElementType,mark_temp,iDummy,iDummy, &
      			   points_temp(1),points_temp(2)

           if(nElementType==1) then
			iTemp1=iTemp1+1
	     elseif(nElementType==2) then
	        goto 100
	     end if
	 end do

100	 nBoundaryEdges=iTemp1
	 nFaces=nElements - nBoundaryEdges

	 rewind 1

     !read out 4 line of comments
     read(1,*) string
	 read(1,*) string
	 read(1,*) string
	 read(1,*) string
       
	 read(1,*) nNodes	

       do i = 1, nNodes
           read(1, *) iDummy,pcoor(i,1),pcoor(i,2),pcoor(i,3)
       enddo
	
	!read out 2 line of comments 
       read(1,*) string
       read(1,*) string
	 
	 read(1,*) nElements

	 do i = 1, nBoundaryEdges 
 		 read(1,*) iDummy,nElementType,iDummy,mark_temp,iDummy,iDummy, &
     			   points_temp(1),points_temp(2)

	     	boundaryEdgeMarkers(i)=mark_temp
			boundaryEdgePoints(i,1)=points_temp(1)
			boundaryEdgePoints(i,2)=points_temp(2)
	 enddo

	 do i = 1, nFaces 
 		 read(1,*) iDummy,nElementType,iDummy, mark_temp,iDummy,iDummy, &
     			   points_temp(1),points_temp(2),points_temp(3)

			facePoints(i,1)=points_temp(1)
			facePoints(i,2)=points_temp(2)
			facePoints(i,3)=points_temp(3)

			faceEdgesNum(i)=3 !GMSH only has triangle elements
			facePoints(i,4)=points_temp(1) !make a loop
 	 enddo

     close(1)

	end subroutine
