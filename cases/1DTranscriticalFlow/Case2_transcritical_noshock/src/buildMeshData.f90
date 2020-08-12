    ! build mesh's topological data which will not change: 
	!                  e.g. face's edges 
    !                       edge's faces
    !                       face area
    !                       face normals
    !                       point's faces

    !edge's points(edgePoints)
    !       in edgePoints, always small number of points in the begining
    !loop over each face to collect not previously stored edges

	subroutine buildMeshData()
	USE COMMON_MODULE
	implicit none
	
	integer*4 i,j,k
    real*8 v1(3), v2(3), v3(3), v4(3), v5(3)
    real*8 w1(2), w2(2), c1(3), c2(3), edgeCenter(3)
    !curEdgeNormal: edge's length normal(norm equals edge's length)
    real*8 curEdgeNormal(3)

    real*8 dhdx, dhdy, dhdz
    integer*4 iIndicator
    integer*4 face1, face2
    integer*4 pointsTemp(ELEDGES+1)
    real*8 dtemp1, dtemp2, temp

	logical bIndicator, bFound

!   Functions out side of main program file
    real*8 dist2points        

	write(*,*) 'Building mesh data ...'

	!setup the elevation of each nodes
	call interp_Elevation

  	  curEdge = 0
      do i = 1, nFaces
         !current face's points
		 do j=1, faceEdgesNum(i)
			pointsTemp(j) = facePoints(i,j)
		 end do
		 pointsTemp(faceEdgesNum(i)+1)=facePoints(i,1)

         if(i.eq.1) then !first face's edges
			do j=1,faceEdgesNum(i)
               curEdge = curEdge + 1
               edgePoints(curEdge,1) = pointsTemp(j)
               edgePoints(curEdge,2) = pointsTemp(j+1)
			   faceEdges(i,j)=curEdge
			end do
         else            !other face's edges
			do j=1,faceEdgesNum(i)
				 !check if the edge is already in the edge list
				 bIndicator = .FALSE.
				 do k = 1, curEdge 
					if((edgePoints(k,1).eq.pointsTemp(j).and.edgePoints(k,2).eq.pointsTemp(j+1)).or.&
					   (edgePoints(k,1).eq.pointsTemp(j+1).and.edgePoints(k,2).eq.pointsTemp(j))) then !already in the list
					   bIndicator = .TRUE.
					   faceEdges(i,j)=k
					   exit
					endif
				 enddo
           
				 if(.not.bIndicator) then !a new edge
					 curEdge = curEdge + 1
					 edgePoints(curEdge,1) = pointsTemp(j)
					 edgePoints(curEdge,2) = pointsTemp(j+1)
					 faceEdges(i,j)=curEdge
				 endif
			end do
         end if
      enddo
      !total number of edges
      nEdges = curEdge

      !for debug
!      write(*,*) 'edgePoints:'
      do i = 1, nEdges
!         write(*,*) edgePoints(i,1), edgePoints(i,2)
      end do

      !edge's faces
      do i = 1, nEdges
         edgeFaces(i, 1) = -1
         edgeFaces(i, 2) = -1

         iIndicator = 0
         !loop over each face to find the owner and neighbour
         do j = 1, nFaces
			do k = 1, faceEdgesNum(j)
				if(faceEdges(j,k).eq.i) then
					iIndicator = iIndicator + 1
					edgeFaces(i,iIndicator) = j
				end if
			end do
		 end do
          
		 if((iIndicator.ge.3).or.(iIndicator.le.0)) then !something bad happended
             write(*,*) 'The edge is shared by more than 2 face or no face!'
		     
             stop
         end if
                 
      enddo

	!build the information of face's edges markers
	binfo = 4     !default is 4-->internal
	pointMarkers = 0 !default is outside of the domain
	do i=1,nFaces
		do j=1,faceEdgesNum(i)
			call checkEdgeMarkers(faceEdges(i,j),binfo(i,j))
			edgeMarkers(faceEdges(i,j)) = binfo(i,j)
			pointMarkers(facePoints(i,j))=1
		end do
	enddo

	!build the information of face's neighbors
	iIndicator=0 !ghost cells counter
	do i=1,nFaces
		do j=1,faceEdgesNum(i)
			!get face's j direction edge first
			curEdge = faceEdges(i,j)

			if(edgeMarkers(curEdge)==1) then      !inlet
				iIndicator=iIndicator+1
				faceNeighbors(i,j) = -1
				boundaryEdgeGhostCells(curEdge) = iIndicator
				ghostCellsEdge(iIndicator)=curEdge
				ghostCellsNeighbor(iIndicator)=i
			elseif(edgeMarkers(curEdge)==2) then  !outlet
				iIndicator=iIndicator+1
				faceNeighbors(i,j) = -2
				boundaryEdgeGhostCells(curEdge) = iIndicator
				ghostCellsEdge(iIndicator)=curEdge
				ghostCellsNeighbor(iIndicator)=i
			elseif(edgeMarkers(curEdge)==3) then  !wall
				iIndicator=iIndicator+1
				faceNeighbors(i,j) = -3
				boundaryEdgeGhostCells(curEdge) = iIndicator
				ghostCellsEdge(iIndicator)=curEdge
				ghostCellsNeighbor(iIndicator)=i
			else                                  !internal
				if(edgeFaces(curEdge,1)/=i) then
					faceNeighbors(i,j)=edgeFaces(curEdge,1)
				else
					faceNeighbors(i,j)=edgeFaces(curEdge,2)
				endif
			endif
		enddo
	enddo

	!Test the total number of ghost cells and boundary edges
	if(iIndicator.ne.nBoundaryEdges)then
		write(*,*) "Ghost cells problem!!"
		
		stop
	endif

	!build the information of inlet faces
	do i=1,nFaces
		!default value is 4(non-inlet)
		pinfo(i) = 4
		
		do j=1,faceEdgesNum(i)
			if(binfo(i,j).eq.1) then
				pinfo(i) = 1
			endif
		end do

	enddo
      
      !for debug
!      write(*,*) 'face Edges:'
      do i = 1, nFaces
!         write(*,*) faceEdges(i,1), faceEdges(i,2), faceEdges(i,3)     
      end do 

      !point's faces
      !for each node, loop over all facePoints
      do i =1, nNodes
         do j = 1, maxnodefaces_
            pointFaces(i,j) = -1
         end do

         iIndicator = 0
         do j = 1, nFaces
			do k = 1, faceEdgesNum(j)
				if(facePoints(j,k).eq.i)then
					iIndicator = iIndicator + 1

					if(iIndicator.gt.maxnodefaces_) then
                     write(*,*) 'node face number greater than  &
                         maxnodefaces_! Increase maxnodefaces_.'
					 
                     stop
                    end if
                  
					pointFaces(i,iIndicator)=j
				end if	
			end do
         end do
		 pointNFaces(i)=iIndicator
      end do      

      !for debug
!      write(*,*) 'node faces'
      do i = 1, nNodes
!         write(*,*) pointFaces(i,1), pointFaces(i,2), pointFaces(i,3),
!     *              pointFaces(i,4), pointFaces(i,5), pointFaces(i,6)
      end do

	  !edge's length (projected on z=0 plane, 2D length)
	  do i = 1, nEdges
		 v1(1) = pcoor(edgePoints(i,1),1)
		 v1(2) = pcoor(edgePoints(i,1),2)
		 v1(3) = 0.0

		 v2(1) = pcoor(edgePoints(i,2),1)
		 v2(2) = pcoor(edgePoints(i,2),2)
		 v2(3) = 0.0
		  
		 edgeLength(i) = dist2points(v1,v2)
	  enddo

      !face area projected on z=0 plane, 2D area
      do i = 1, nFaces
		  face2DArea(i)=0.0
		  do j=1,faceEdgesNum(i)
			 face2DArea(i)=face2DArea(i)+pcoor(facePoints(i,j),1)*pcoor(facePoints(i,j+1),2)
			 face2DArea(i)=face2DArea(i)-pcoor(facePoints(i,j),2)*pcoor(facePoints(i,j+1),1)
		  end do
		
		  face2DArea(i)=face2DArea(i)/2.0

          if(face2DArea(i).lt.0.0) then
              write(*,*) 'Negative 2D area for face:', i
              face2DArea(i) = -1.0*face2DArea(i)
          end if
      enddo

	  !face's edges normal vector (pointing out of the cell)
	  do i=1, nFaces
	     !face center
		 faceCenters(i,1)=0.0
		 faceCenters(i,2)=0.0
		 faceCenters(i,3)=0.0
		 do j=1,faceEdgesNum(i)
			faceCenters(i,1)=faceCenters(i,1)+pcoor(facePoints(i,j),1)
			faceCenters(i,2)=faceCenters(i,2)+pcoor(facePoints(i,j),2)
			faceCenters(i,3)=faceCenters(i,3)+pcoor(facePoints(i,j),3)
		 end do
		 faceCenters(i,1)=faceCenters(i,1)/faceEdgesNum(i)
		 faceCenters(i,2)=faceCenters(i,2)/faceEdgesNum(i)
		 faceCenters(i,3)=faceCenters(i,3)/faceEdgesNum(i)

		 c1(1)=faceCenters(i,1)	
		 c1(2)=faceCenters(i,2)
		 c1(3)=faceCenters(i,3)

		do j=1,faceEdgesNum(i)   !for each edge
			curEdge=faceEdges(i,j)
            !edge vector: from start to end
            v1(1) = pcoor(edgePoints(curEdge, 1), 1) -  &
                    pcoor(edgePoints(curEdge, 2), 1)
            v1(2) = pcoor(edgePoints(curEdge, 1), 2) -  &
                    pcoor(edgePoints(curEdge, 2), 2)
            v1(3) = pcoor(edgePoints(curEdge, 1), 3) -  &
                    pcoor(edgePoints(curEdge, 2), 3)

            !edge center
            edgeCenterCoor(curEdge,1) = 0.5*(pcoor(edgePoints(curEdge, 1), 1) +  &
                         pcoor(edgePoints(curEdge, 2), 1) )
            edgeCenterCoor(curEdge,2) = 0.5*(pcoor(edgePoints(curEdge, 1), 2) +  &
                         pcoor(edgePoints(curEdge, 2), 2) )
			edgeCenterCoor(curEdge,3) = 0.5*(pcoor(edgePoints(curEdge, 1), 3) +  &
                         pcoor(edgePoints(curEdge, 2), 3) )
        
            !aux vector: from face center to edge center (in 2D)
            v3(1) = edgeCenterCoor(curEdge,1) - c1(1)
            v3(2) = edgeCenterCoor(curEdge,2) - c1(2)
            v3(3) = 0

            !current edge's normal, the norm is the length of the edge
            curEdgeNormal(1) = -v1(2)
            curEdgeNormal(2) = v1(1)
            curEdgeNormal(3) = 0

			!normalize the vector
			call vec_normalize(curEdgeNormal)

            !adjust the normal vector to point out of the triangle
            call vec_dot(curEdgeNormal, v3, temp)
            faceEdgeNormals(i,j,1) = sign(1.0D0, temp)*curEdgeNormal(1)
            faceEdgeNormals(i,j,2) = sign(1.0D0, temp)*curEdgeNormal(2)
            faceEdgeNormals(i,j,3) = sign(1.0D0, temp)*curEdgeNormal(3)
		enddo
	  enddo

	  !calculate face's nearest distance from center to vertices (2D)
	  do i=1, nFaces
		faceMinR(i)=VLARGE

        v1(1) = faceCenters(i,1)
        v1(2) = faceCenters(i,2)
        v1(3) = 0.0D0
		
		do j=1, faceEdgesNum(i)
			v2(1) = pcoor(facePoints(i,j),1)
			v2(2) = pcoor(facePoints(i,j),1)
			v2(3) = 0.0D0
		enddo

		faceMinR(i)=min(faceMinR(i),dist2points(v1,v2))

	  enddo

  	  !calculate faceSubArea
	  do i = 1, nFaces
        v1(1) = faceCenters(i,1)
        v1(2) = faceCenters(i,2)
        v1(3) = 0

		do j = 1, faceEdgesNum(i)
          v2(1) = pcoor(edgePoints(faceEdges(i, j),1), 1)
          v2(2) = pcoor(edgePoints(faceEdges(i, j),1), 2)
          v2(3) = 0

          v3(1) = pcoor(edgePoints(faceEdges(i, j),2), 1)
          v3(2) = pcoor(edgePoints(faceEdges(i, j),2), 2)
          v3(3) = 0

          !signed area: countclock positive
          !             clockwise  negative
          faceSubArea(i,j) = 0.5*(-v2(1)*v1(2)+v3(1)*v1(2)+v1(1)*v2(2)  &
                                  -v3(1)*v2(2)-v1(1)*v3(2)+v2(1)*v3(2))

          if(faceSubArea(i,j).lt.0.0) then
              !write(*,*) 'Negative 2D area for face:', i
              faceSubArea(i,j) = -1.0*faceSubArea(i,j)
          end if
		end do
	  end do

	  !calculate the mesh information which will change with the scour process
	  call updateMeshData
	end subroutine


!	calculate the mesh information which will change with the updating of bed
!   elevation
	subroutine updateMeshData()
	USE COMMON_MODULE
	implicit none

	integer i,j
	real*8 curEdgeNormal(3),temp


	!update bed slope
	Sox=0.0
	Soy=0.0

	do i = 1, nFaces
	   do j = 1, faceEdgesNum(i)	    
		temp=(pcoor(edgePoints(faceEdges(i,j),1),3)+pcoor(edgePoints(faceEdges(i,j),2),3))/2* &
		     faceEdgeNormals(i,j,1)*edgeLength(faceEdges(i,j))
		Sox(i)=Sox(i)-1.0/face2DArea(i)*temp
	
		temp=(pcoor(edgePoints(faceEdges(i,j),1),3)+pcoor(edgePoints(faceEdges(i,j),2),3))/2* &
		     faceEdgeNormals(i,j,2)*edgeLength(faceEdges(i,j))
		Soy(i)=Soy(i)-1.0/face2DArea(i)*temp
	   end do
	end do

	!update edge center coordinates
	do i = 1, nEdges
       edgeCenterCoor(i,1) = 0.5*(pcoor(edgePoints(i, 1), 1) +  &
		                          pcoor(edgePoints(i, 2), 1) )
       edgeCenterCoor(i,2) = 0.5*(pcoor(edgePoints(i, 1), 2) +  &
                                  pcoor(edgePoints(i, 2), 2) )
       edgeCenterCoor(i,3) = 0.5*(pcoor(edgePoints(i, 1), 3) +  &
                                  pcoor(edgePoints(i, 2), 3) )
	end do

	return
	end

!   get the edge marker for a given edge 
	subroutine checkEdgeMarkers(edge, boundaryMarker)
	USE COMMON_MODULE
	implicit none

	integer*4 edge, boundaryMarker

	integer*4 p1,p2
	integer*4 i
	logical bFound

	bFound=.false.

	p1=edgePoints(edge,1)
	p2=edgePoints(edge,2)

	!first, look at the boundary edge set
	do i=1, nBoundaryEdges
		if(((p1.eq.boundaryEdgePoints(i,1)).and.(p2.eq.boundaryEdgePoints(i,2))).or.&
		   ((p1.eq.boundaryEdgePoints(i,2)).and.(p2.eq.boundaryEdgePoints(i,1)))) then
			bFound=.true.
			boundaryMarker=boundaryEdgeMarkers(i)
			exit
		end if
	end do

	!if not in the boundary edge set, then it's an internal (code 4)
	if(.not.bFound) then
		boundaryMarker=4
	end if

	end subroutine


!   interpolate the bed elevation according to DEM file
	subroutine interp_Elevation()
	USE COMMON_MODULE,only:pcoor,nNodes
	implicit none

	integer i
	real*8 r2, r0, r

	!set up the bottom elevation for the circular basin test case
!	r0=192
!	do i=1,nNodes
!		r=dsqrt((pcoor(i,1)-192)**2+(pcoor(i,2)-192)**2)
!		pcoor(i,3)=-(0.5+dsqrt(dabs(0.5-0.5*r/r0)))/1.3
!	enddo 

	!set up the bottom elevation for the long channel case
!	do i=1,nNodes
!		if(pcoor(i,1).le.700.and.pcoor(i,1).ge.500.and. &
!		   pcoor(i,2).le.600.and.pcoor(i,2).ge.400) then
!		   pcoor(i,3)=dsin(3.1415926*(pcoor(i,1)-500)/200)**2*&
!					  dsin(3.1415926*(pcoor(i,2)-400)/200)**2
!		end if
!	enddo 

    !set up for the transcritical flow case
	do i=1,nNodes
	   pcoor(i,3)=max(0.0D0, 0.2-0.05*(pcoor(i,1)-10)**2)
	enddo 

	!rect basin with a hump case
!	do i=1,nNodes
!		r2=-50*(pcoor(i,1)-0.5)**2-50*(pcoor(i,2)-0.5)**2
!		pcoor(i,3)=0.5*dexp(r2)
!	enddo 


	!1D Gaussian hump case
!	do i=1,nNodes
!		pcoor(i,3)=-6+2*dexp(-0.01*(pcoor(i,1)-150)**2)
!	end do


	!set up the bottom elevation for the 90 degree bend case
!	do i=1,nNodes
!		if(pcoor(i,1).ge.2.39) then
!		   pcoor(i,3)=0.33
!		end if
!	enddo 

	end subroutine
