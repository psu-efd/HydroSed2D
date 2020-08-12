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

!   sediment transport: one step
	subroutine sediment()
	USE COMMON_MODULE
	implicit none

	!calculate shear stress at faces
	call calc_shear

	!calculate steepest slope
	call calc_steepestSlope

    !calculate critical Shields number
	call calc_criticalShieldsNumber

	!calculate face Shields number
	call calc_ShieldsNumber

	!calculate sediment flux at face center and edge center
	call calc_sedimentFlux

    !calculate divergence of sediment flux at face center
	call calc_divSedimentFlux

    !calculate face center elevation change according to Exner equation
	call calc_faceElevChange 

	!update the nodal value of elevation 
	call calc_nodeNewElevation

	!update the mesh information
	call updateMeshData

	return
	end

!   read sediment properties
    subroutine readSedimentProperties()
	USE COMMON_MODULE
	implicit none

	character*100 string

	if(platform.eq.1) then  !windows platform
		tempfilename=trim(trim(dirname) // '\input\sediment.dat')
	else                    !linux platform 
	    tempfilename=trim(trim(dirname) // '/input/sediment.dat')
	endif

    open(1, file = tempfilename, status = 'old', iostat = ios)
    if ( ios /= 0 ) then
           ierror = 1
           write ( *, '(a)' ) ' '
           write ( *, '(a)' ) 'Could not open the input file: sediment.dat'
           stop
    end if

    read(1, *) rhos, string
    read(1, *) rhow, string 
    read(1, *) diam, string 
    read(1, *) porosity, string 
    read(1, *) thita_cri_b, string 
    read(1, *) thita_cri_s, string 
    read(1, *) angleofRepose, string
    read(1, *) C_slope, string 
    read(1, *) sigmac, string
	read(1, *) sedTransRate, string
	read(1, *) GrassConst(1),GrassConst(2),string
       
    close(1)

	!submerged relative density
	Rs=(rhos-rhow)/rhow

	return
	end 


!	calculate shear stress at nodes using Manning-Strickler formula
	subroutine calc_shear()
	USE COMMON_MODULE
	implicit none

	integer i
	real*8 sfx,sfy
	
	!face shear stress
	do i = 1, nFaces
		sfx=zaolv**2*UM(i)*dsqrt(UM(i)**2+VN(i)**2)/Q1(i)**(4.0/3.0)
		sfy=zaolv**2*VN(i)*dsqrt(UM(i)**2+VN(i)**2)/Q1(i)**(4.0/3.0)		                  
		faceShearStress(i,1)=rhow*g*Q1(i)*sfx
		faceShearStress(i,2)=rhow*g*Q1(i)*sfy
    end do

	return
	end

!	calculate face Shields number
	subroutine calc_ShieldsNumber()
	USE COMMON_MODULE
	implicit none

	integer i
	real*8 temp

	do i = 1, nFaces
		temp=dsqrt(faceShearStress(i,1)**2+faceShearStress(i,2)**2)
		faceShieldsNumber(i)= temp/(rhow*g*Rs*diam)
	end do
	

	return
	end 


!	calculate face critical Shields number (slope effect)
	subroutine calc_criticalShieldsNumber()
	USE COMMON_MODULE
	implicit none

	integer i
	real*8 temp,v1(3),v2(3),flowdir

	do i = 1, nFaces
       ! sign function to judge the flow direction: uphill or downhill
       v1(1)=faceShearStress(i,1)
       v1(2)=faceShearStress(i,2)
       v1(3)=0.0D0
       v2(1)=faceSteepestSlope(i,1)
       v2(2)=faceSteepestSlope(i,2)
       v2(3)=faceSteepestSlope(i,3)

       call vec_dot(v1, v2, flowDir)
           
       faceCriticalShieldsNumber(i) = thita_cri_b* &
             dsin(angleofRepose+dsign(1.0D0,flowDir)*faceSlopeAngle(i))/ &
                              dsin(angleofRepose)
	end do
	
	return
	end 

!   calculate sediment flux
	subroutine calc_sedimentFlux()
	USE COMMON_MODULE
	implicit none

	integer i,face1,face2
	real*8 dzdx,dzdy,temp,dtemp1,dtemp2
	real*8 v1(2),c1(3),c2(3),edgeCenter(3)
	real*8 dist2points
	real*8 utemp,vtemp

    !calculate sediment flux at face centers
    do i =1, nFaces
	  if(sedTransRate.eq.1) then        !MPM formula
         if(faceShieldsNumber(i).gt.faceCriticalShieldsNumber(i)) then
	         !flux at flat bed: qflat
			 qflat = 8*dsqrt(g*Rs*diam**3)*(faceShieldsNumber(i)-faceCriticalShieldsNumber(i))**(1.5)

			 !slope effect of the sediment flux
			 dzdx = -Sox(i) 
			 dzdy = -Soy(i)

			 v1(1) = faceShearStress(i, 1)
			 v1(2) = faceShearStress(i, 2)
			 temp=dsqrt(faceShearStress(i,1)**2+faceShearStress(i,2)**2)

			 faceBedLoadFlux(i,1)=qflat*faceShearStress(i,1)/(temp+VSMALL)-C_slope*dabs(qflat)*dzdx
			 faceBedLoadFlux(i,2)=qflat*faceShearStress(i,2)/(temp+VSMALL)-C_slope*dabs(qflat)*dzdy
		  else
 			 faceBedLoadFlux(i,1)=0.0
			 faceBedLoadFlux(i,2)=0.0
		  end if
	   else if(sedTransRate.eq.2) then   !Grass formula
 			 faceBedLoadFlux(i,1)=GrassConst(1)*UM(i)*(UM(i)**2+VN(i)**2)**(0.5*(GrassConst(2)-1))
			 faceBedLoadFlux(i,2)=GrassConst(1)*VN(i)*(UM(i)**2+VN(i)**2)**(0.5*(GrassConst(2)-1))

!			 utemp=10/(0-faceCenters(i,3))
!			 vtemp=0.0
! 			 faceBedLoadFlux(i,1)=GrassConst(1)*utemp*(utemp**2+vtemp**2)**(0.5*(GrassConst(2)-1))
!			 faceBedLoadFlux(i,2)=GrassConst(1)*vtemp*(utemp**2+vtemp**2)**(0.5*(GrassConst(2)-1))
	   else
			stop 'Sediment transport formula not supported!'
	   end if
     enddo
     
     !interpolate the face center values to edges
     do i = 1, nEdges
         if(edgeFaces(i,2).eq.-1) then !boundary edges
            edgeBedLoadFlux(i,1) = faceBedLoadFlux(edgeFaces(i,1),1)
            edgeBedLoadFlux(i,2) = faceBedLoadFlux(edgeFaces(i,1),2)
         else                          !internal edges
            !face1: owner   face2: neighbour
            face1 = edgeFaces(i,1)
            face2 = edgeFaces(i,2)
            !c1: owner's center   c2: neighbour's center
            c1(1) = faceCenters(face1,1)
            c1(2) = faceCenters(face1,2)
            c1(3) = 0.0D0 

            c2(1) = faceCenters(face2,1)
            c2(2) = faceCenters(face2,2)
            c2(3) = 0.0D0 

            !current edge center
            edgeCenter(1) = edgeCenterCoor(i,1)
            edgeCenter(2) = edgeCenterCoor(i,2)
            edgeCenter(3) = 0.0D0

            !dtemp1:distance between owner's center to edge center
            !dtemp2:distance between neighbour's center to edge center
            dtemp1 = dist2points(c1, edgeCenter)
            dtemp2 = dist2points(c2, edgeCenter)
           
            !weighted average by inverse distance
            edgeBedLoadFlux(i, 1) = (faceBedLoadFlux(face1,1)*dtemp2+&
                                     faceBedLoadFlux(face2,1)*dtemp1)/&
                                    (dtemp1+dtemp2)
            edgeBedLoadFlux(i, 2) = (faceBedLoadFlux(face1,2)*dtemp2+&
                                     faceBedLoadFlux(face2,2)*dtemp1)/&
                                    (dtemp1+dtemp2)
         end if
     enddo


	return 
	end


!	calculate divergence of sediment flux at face centers
	subroutine calc_divSedimentFlux()
	USE COMMON_MODULE
	implicit none

	integer i,j,pos
    real*8 c1(3),v1(3),v2(3),v3(3),v4(3),curEdgeNormal(3)
	real*8 temp
	integer edgePositionInFace

    do i = 1, nFaces
         divFaceBedLoadFlux(i) = 0.0

         ! for each edge
         do j = 1, faceEdgesNum(i)
            curEdge = faceEdges(i, j)

			pos=edgePositionInFace(i,curEdge)
			if(pos==-1)then
				write(*,*) 'Something wrong!'
				stop
			else
				curEdgeNormal(1) = faceEdgeNormals(i,pos,1)
				curEdgeNormal(2) = faceEdgeNormals(i,pos,2)
				curEdgeNormal(3) = 0.0				
			end if

            !current edge sediment flux vector(expand to 3D vector for vec_dot)
            v4(1) = edgeBedLoadFlux(curEdge, 1)  !x
            v4(2) = edgeBedLoadFlux(curEdge, 2)  !y 
            v4(3) = 0.0                          !z

            !contribution of divergence from current edge integration
            call vec_dot(v4, curEdgeNormal, temp)
            divFaceBedLoadFlux(i) = divFaceBedLoadFlux(i) + temp*edgeLength(curEdge)      
         end do      
	  end do

	return
	end

    
!   calculate the face elevation change according to Exner equation
!   and updating the faceCenters elevation
	subroutine calc_faceElevChange()
	USE COMMON_MODULE
	implicit none

	integer i

    do i = 1, nFaces
       faceElevationChange(i) = -sedInterval*dt/(1-porosity)*divFaceBedLoadFlux(i)/face2DArea(i)
	   faceCenters(i,3) = faceCenters(i,3) + faceElevationChange(i)
    end do

	return
	end

!   calculate node new elevations
	subroutine calc_nodeNewElevation()
	USE COMMON_MODULE
	implicit none

	integer i,j

	!interpolate the face elevation change to node values
	call cellCenterToNodes(faceElevationChange,nodeElevationChange)

	do i = 1, nNodes
		pcoor(i,3) = pcoor(i,3) + nodeElevationChange(i)
	end do

	!recalculat the face center elevation to syncronize after node coor change
    do i=1, nFaces
	     !face center
		 faceCenters(i,3)=0.0
		 do j=1,faceEdgesNum(i)
			faceCenters(i,3)=faceCenters(i,3)+pcoor(facePoints(i,j),3)
		 end do
		 faceCenters(i,3)=faceCenters(i,3)/faceEdgesNum(i)
	end do

	return
	end

!   calculate steepest slope
	subroutine calc_steepestSlope()
	USE COMMON_MODULE, ONLY:Sox,Soy,faceSlopeAngle,faceSteepestSlope,nFaces
	implicit none

	integer i
	real*8 A, B, D, norm

    !calculate steepest slope and slope angle given unit normal of the face
    do i = 1, nFaces
          A = -Sox(i)
          B = -Soy(i)
          D = 1.0

		  norm = dsqrt(A**2+B**2+D**2)
		  !normalize
		  A = A/norm
		  B = B/norm
		  D = D/norm
      
          faceSlopeAngle(i) = dasin(A**2+B**2)
          faceSteepestSlope(i, 1) = A*D/dsqrt(A**2+B**2+1E-16)
          faceSteepestSlope(i, 2) = D*B/dsqrt(A**2+B**2+1E-16)
          faceSteepestSlope(i, 3) = -dsqrt(A**2+B**2)
    end do


	return
	end
