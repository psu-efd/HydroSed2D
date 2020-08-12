      ! build mesh data: e.g. face's edges 
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

	logical bIndicator, bFound
	
      
	curEdge = 0
      do i = 1, nFaces
         !current face's three points
         curPoint1 = facePoints(i,1)
         curPoint2 = facePoints(i,2)
         curPoint3 = facePoints(i,3)
         !order those in ascending order
         call orderThreeInt(curPoint1, curPoint2, curPoint3)

         if(i.eq.1) then !first face's three edges
             curEdge = curEdge + 1
             edgePoints(curEdge,1) = curPoint1
             edgePoints(curEdge,2) = curPoint2
             curEdge = curEdge + 1
             edgePoints(curEdge,1) = curPoint1
             edgePoints(curEdge,2) = curPoint3
             curEdge = curEdge + 1
             edgePoints(curEdge,1) = curPoint2
             edgePoints(curEdge,2) = curPoint3
         else            !other face's three edges
             !check if first edge is already in the edge list
             bIndicator = .FALSE.
             do k = 1, curEdge 
                if(edgePoints(k,1).eq.curPoint1.and.
     *             edgePoints(k,2).eq.curPoint2) then !already in the list
                   bIndicator = .TRUE.
                   exit
                endif
             enddo
           
             if(.not.bIndicator) then !a new edge
                 curEdge = curEdge + 1
                 edgePoints(curEdge,1) = curPoint1
                 edgePoints(curEdge,2) = curPoint2
             endif
     
             !check if second edge is already in the edge list
             bIndicator = .FALSE.
             do k = 1, curEdge
                if(edgePoints(k,1).eq.curPoint1.and.
     *             edgePoints(k,2).eq.curPoint3) then !already in the list
                   bIndicator = .TRUE.
                   exit
                endif
             enddo

             if(.not.bIndicator) then !a new edge
                 curEdge = curEdge + 1
                 edgePoints(curEdge,1) = curPoint1
                 edgePoints(curEdge,2) = curPoint3
             endif

             !check if third edge is already in the edge list
             bIndicator = .FALSE.
             do k = 1, curEdge
                if(edgePoints(k,1).eq.curPoint2.and.
     *             edgePoints(k,2).eq.curPoint3) then !already in the list
                   bIndicator = .TRUE.
                   exit
                endif
             enddo

             if(.not.bIndicator) then !a new edge
                 curEdge = curEdge + 1
                 edgePoints(curEdge,1) = curPoint2
                 edgePoints(curEdge,2) = curPoint3
             endif
         end if
      enddo
      !total number of edges
      nEdges = curEdge

      !for debug
      write(*,*) 'edgePoints:'
      do i = 1, nEdges
!         write(*,*) edgePoints(i,1), edgePoints(i,2)
      end do

      !face's edges
      do i = 1, nFaces
         faceEdges(i,1) = -1
         faceEdges(i,2) = -1
         faceEdges(i,3) = -1

         !current face's three points
         curPoint1 = facePoints(i,1)
         curPoint2 = facePoints(i,2)
         curPoint3 = facePoints(i,3)
         !order those in ascending order
         call orderThreeInt(curPoint1, curPoint2, curPoint3)

         !find first edge number
         bFound = .FALSE.
         do j = 1, nEdges
 !           write(*,*) edgePoints(j,1),edgePoints(j,2)

            if(edgePoints(j,1).eq.curPoint1.and.
     *         edgePoints(j,2).eq.curPoint2) then !got you
               faceEdges(i,1) = j
               bFound = .TRUE.
               exit
            end if
         enddo
         if(.not.bFound) then
            write(*,*) 'Did not find the edge number for face:', i
            stop
         end if

         !find second edge number
         bFound = .FALSE.
         do j = 1, nEdges
            if(edgePoints(j,1).eq.curPoint1.and.
     *         edgePoints(j,2).eq.curPoint3) then !got you
               faceEdges(i,2) = j
               bFound = .TRUE.
               exit
            end if
         enddo
         if(.not.bFound) then
            write(*,*) 'Did not find the edge number for face:', i
            stop
         end if

         !find third edge number
         bFound = .FALSE.
         do j = 1, nEdges
            if(edgePoints(j,1).eq.curPoint2.and.
     *         edgePoints(j,2).eq.curPoint3) then !got you
               faceEdges(i,3) = j
               bFound = .TRUE.
               exit
            end if
         enddo
         if(.not.bFound) then
            write(*,*) 'Did not find the edge number for face:', i
            stop
         end if
      enddo
      
      !for debug
      write(*,*) 'face Edges:'
      do i = 1, nFaces
!         write(*,*) faceEdges(i,1), faceEdges(i,2), faceEdges(i,3)     
      end do 

      !edge's faces
      do i = 1, nEdges
         edgeFaces(i, 1) = -1
         edgeFaces(i, 2) = -1

         iIndicator = 0
         !loop over each face to find the owner and neighbour
         do j = 1, nFaces
            if(faceEdges(j, 1).eq.i.or.
     *         faceEdges(j, 2).eq.i.or.
     *         faceEdges(j, 3).eq.i) then
               iIndicator = iIndicator +1
               if(iIndicator.ge.3) then !something bad happended
                                        !the edge is shared by more than 2 faces
                  write(*,*) 'The edge is shared by more than 2 face!'
                  stop
               end if
               
               edgeFaces(i,iIndicator) = j

            end if
         enddo
      enddo

      !for debug
      write(*,*) 'edge Faces:'
      do i = 1, nEdges
!         write(*,*) edgeFaces(i,1), edgeFaces(i,2)
      end do

      !face area projected on y=0 plane, 2D area
      !should be changed to z=0 plane in femco
      do i = 1, nFaces
          v1(1) = pcoor(facePoints(i, 1), 1)
          v1(2) = 0.0
          v1(3) = pcoor(facePoints(i, 1), 3)

          v2(1) = pcoor(facePoints(i, 2), 1)
          v2(2) = 0.0
          v2(3) = pcoor(facePoints(i, 2), 3)

          v3(1) = pcoor(facePoints(i, 3), 1)
          v3(2) = 0.0
          v3(3) = pcoor(facePoints(i, 3), 3)

          !signed area: countclock positive
          !             clockwise  negative
          face2DArea(i) = 0.5*(-v2(1)*v1(3)+v3(1)*v1(3)+v1(1)*v2(3)
     *                         -v3(1)*v2(3)-v1(1)*v3(3)+v2(1)*v3(3))

          if(face2DArea(i).lt.0.0) then
!              write(*,*) 'Negative 2D area for face:', i
              face2DArea(i) = -1.0*face2DArea(i)
          end if
      enddo

      !face's unit normal vector
      !not sure about the direction!!!!
      do i = 1, nFaces
         !v1: edges vector 1->2
         v1(1) = pcoor(facePoints(i,1),1) - pcoor(facePoints(i,2),1)
         v1(2) = pcoor(facePoints(i,1),2) - pcoor(facePoints(i,2),2)
         v1(3) = pcoor(facePoints(i,1),3) - pcoor(facePoints(i,2),3)
         !v2: edges vector 1->3
         v2(1) = pcoor(facePoints(i,1),1) - pcoor(facePoints(i,3),1)
         v2(2) = pcoor(facePoints(i,1),2) - pcoor(facePoints(i,3),2)
         v2(3) = pcoor(facePoints(i,1),3) - pcoor(facePoints(i,3),3)

         !take cross pruduction of v1 X v2 --> v3
         call vec_cross(v1, v2, v3)
         !normalize v3
         call vec_normalize(v3)

         faceUnitNormals(i,1) = v3(1)
         faceUnitNormals(i,2) = v3(2)
         faceUnitNormals(i,3) = v3(3)
      enddo
 
      !for debug
      write(*,*) 'faceUnitNormals'
      do i = 1, nFaces
!         write(*,*) faceUnitNormals(i,1), faceUnitNormals(i,2),
!     *              faceUnitNormals(i,3)
      end do

      !point's faces
      !for each node, loop over all facePoints
      do i =1, nNodes
         do j = 1, maxnodefaces_
            pointFaces(i,j) = -1
         end do

         iIndicator = 0
         do j = 1, nFaces
            if(facePoints(j,1).eq.i.or.
     *         facePoints(j,2).eq.i.or.
     *         facePoints(j,3).eq.i) then
                  iIndicator = iIndicator + 1
                  if(iIndicator.gt.maxnodefaces_) then
                     write(*,*) 'node face number greater than
     *                   maxnodefaces_! Increase maxnodefaces_.'
                     stop
                  end if
                  pointFaces(i,iIndicator)=j
            end if
         end do
      end do      

      !for debug
      write(*,*) 'node faces'
      do i = 1, nNodes
!         write(*,*) pointFaces(i,1), pointFaces(i,2), pointFaces(i,3),
!     *              pointFaces(i,4), pointFaces(i,5), pointFaces(i,6)
      end do
	end subroutine