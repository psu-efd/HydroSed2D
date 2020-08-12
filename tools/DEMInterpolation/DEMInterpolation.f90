!****************************************************************************
!
!  PROGRAM: DEMInterpolation
!
!  PURPOSE: Interpolate the DEM data to the 2D mesh.
!
!****************************************************************************

	program DEMInterpolation
    USE COMMON_MODULE

	implicit none

    integer i,j

    integer iDummy,nElementType
	integer nFileType  !1: GMSH  2: Gambit Neutral

        integer, parameter::nDEMPoints=1433142
	
        real*8 points(nDEMPoints,3),tmp
	CHARACTER*80 DUMMY,string

	real*8 time1, time2

	!define some interpolation constants
	real*8 radius
	integer pointsClose(nDEMPoints) !the DEM point which is closest to current mesh node
    integer upLevel

	CHARACTER*80 NEUFILE,UNVFILE,CONTRINFO,GFILE,HED,PROG,REC5,&
			     ENDOFSECTION,MTHYEAR,GROUP,ELEMENT,MATERIAL,NSFLAGS,&
			     ELMAT,TODAY,MATNAM,boundaryName
	INTEGER NGRPS,NBSETS,NDFCD,NDFVL,ND,NE,&
			NTYPE,NDP,NODE(30),NGP,NELGP,MTYP,NFLAGS,ISOLVE(10),&
			NELT(1000000),REDATE(8),MATNUM,NODET,IType,NEntry,NValues,boundaryType
	integer bFaceTemp, bFaceTypeTemp, bFaceDirTemp,pTemp1,pTemp2
	real*8  boundarFaceTemp 
	REAL(8) REVL

	!information of edge for faces
	integer triInfo(2,3),quadInfo(2,4)
	data triInfo/1,2,2,3,3,1/,quadInfo/1,2,2,3,3,4,4,1/	

    call cpu_time(time1)

	iDummy=1

	!read in the DEM data
	open(2,FILE='dem.dat',STATUS='UNKNOWN')

	do i=1,nDEMPoints
	  read(2,*) points(i,1),points(i,2),points(i,3)
	enddo

	write(*,*) 'Successfully read the DEM data file!'

	close(2)
    
	nFileType=2
	!read the msh file
	if(nFileType.eq.1) then
	    call readGMSHMesh
	else 
	    call readGambitMesh
	endif

    write(*,*) 'Successfully read the mesh file!'

	!interpolate the elevation data to the grid points; use OpenMP
!$omp parallel do private ( i, upLevel, radius, tmp ) shared ( nNodes, pointsClose, nDEMPoints, pcoor, points )
    do i=1, nNodes
	   write(*,*) i, ' out of ', nNodes, ' nodes'
	   upLevel = 0
       radius=1E24
	   pointsClose(i)=-1  !the point close to current node
       do j=1,nDEMPoints
	       upLevel = upLevel + 1
	       tmp=sqrt((points(j,1)-pcoor(i,1))**2+(points(j,2)-pcoor(i,2))**2)
           if(tmp.lt.radius) then !found a closer DEM point
		      radius = tmp
			  pointsClose(i)=j
		   endif
	   enddo

	   if(pointsClose(i).gt.0) then !it has a DEM point close to it
		  pcoor(i,3) = -points(pointsClose(i),3)
	   else
	      write(*,*) 'Node ',i, ' has no depth data.'
	      pcoor(i,3) = -9999
	   endif

	   if(upLevel.gt.18000) then !already went beyond many levels up, quit searching for this node
	      continue
	   endif 
	enddo
!$omp end parallel do

	!write out the depth interpolated data
	open(9,FILE='mesh_depth.dat',STATUS='UNKNOWN')

    write(9,*) 'TITLE = SWE Solver Water Surface Elevation'
	write(9,*) 'VARIABLES = "X", "Y", "Z"'

!	write(9, *) 'ZONE, DATAPACKING=POINT, N=', nNodes, ',E=', nFaces,' ,ZONETYPE=FETRIANGLE'
    write(9, *) 'ZONE, DATAPACKING=POINT, N=', nNodes, ',E=', nFaces,' ,ZONETYPE=FEQUADRILATERAL'
	
	do i=1,nNodes
	  write(9,*) pcoor(i,1),pcoor(i,2),pcoor(i,3)
	enddo

	do i=1,nFaces
	  write(9, *) (facePoints(i,j),j=1,faceEdgesNum(i))
	enddo

	close(9)


    if(nFileType.eq.1) then 
		!write the new GMSH file
		tempfilename='StClaireR_Coarse_new.msh'
		   
		 open(1, file = tempfilename, status = 'unknown', iostat = ios)
		 if ( ios /= 0 ) then
			   ierror = 1
			   write ( *, '(a)' ) ' '
			   write ( *, '(a)' ) 'Could not open the input file'
           
			   stop
		 end if

		 !write out 4 line of comments
		 write(1,*) "TTT"
		 write(1,*) "TTT"
		 write(1,*) "TTT"
		 write(1,*) "TTT"
       
		 write(1,*) nNodes

		 do i = 1, nNodes
			write(1, '(I2, F15.5, F15.5, F12.5)') iDummy,pcoor(i,1),pcoor(i,2),pcoor(i,3)
		 enddo
		
		!write out 2 line of comments 
		 write(1,*) 'TTT'
		 write(1,*) 'TTT'

		!write elements: lines and triangles
		 write(1,*) nElements

		 nElementType=1
		 do i = 1, nBoundaryEdges 
 			 write(1,'(3I2,I2,2I2,2I5)') iDummy,nElementType,iDummy,boundaryEdgeMarkers(i),iDummy,iDummy, &
     				   boundaryEdgePoints(i,1),boundaryEdgePoints(i,2)
		 enddo

		 nElementType=2
		 do i = 1, nFaces 
 			 write(1,'(6I2,3I5)') iDummy,nElementType,iDummy,iDummy,iDummy,iDummy, &
     				   facePoints(i,1),facePoints(i,2),facePoints(i,3)
 		 enddo

		 close(1)
	else
	    !write the new Gambit Neutral file
		tempfilename='StClairRiverUpperFine_new.neu'
		!tempfilename='StClairRiverUpper_new.neu'

		open(9, file = tempfilename, status = 'unknown', iostat = ios)
		if ( ios /= 0 ) then
			 ierror = 1
			 write ( *, '(a)' ) ' '
			 write ( *, '(a)' ) 'Could not open the input file'
         
			 stop
		end if

        !read the original file again, only changed is the vertical coordinate of points
		tempfilename='StClairRiverUpperFine.neu'
		!tempfilename='StClairRiverUpper.neu'

		open(10, file = tempfilename, status = 'old', iostat = ios)
		if ( ios /= 0 ) then
			 ierror = 1
			 write ( *, '(a)' ) ' '
			 write ( *, '(a)' ) 'Could not open the input file'
         
			 stop
		end if


	!	Read Gambit neutral file
	!	Read the header
		READ(10,100) CONTRINFO
		WRITE(9,100) CONTRINFO
		READ(10,100) GFILE
		WRITE(9,100) GFILE
		READ(10,110) HED
		WRITE(9,110) HED
		READ(10,110) PROG
		WRITE(9,110) PROG
		READ(10,110) MTHYEAR
		WRITE(9,110) MTHYEAR
		READ(10,110) REC5
		WRITE(9,110) REC5
		READ(10,130) nNodes,nFaces,NGRPS,NBSETS,NDFCD,NDFVL
		WRITE(9,130) nNodes,nFaces,NGRPS,NBSETS,NDFCD,NDFVL
		READ(10,100) ENDOFSECTION
		WRITE(9,100) ENDOFSECTION
		100 FORMAT(A40)
		110 FORMAT(A80)
		120 FORMAT(A44,F5.2)
		130 FORMAT(6I10)
		140 FORMAT(5X,'nNodes',5X,'nFaces',5X,'NGRPS',4X,'NBSETS',5X,'NDFCD',5X,'NDFVL')
		
	!	Read and write nodal coordinates
		READ(10,'(A30)') CONTRINFO
		WRITE(9,'(A30)') CONTRINFO
		ND=1
		DO WHILE(ND.LT.nNodes)
			READ(10,150) ND,pcoor(ND,1),pcoor(ND,2),tmp  !intentionally don't read pcoor(ND,3)
			WRITE(9,150) ND,pcoor(ND,1),pcoor(ND,2),pcoor(ND,3)
		END DO
		150 FORMAT(I10,1X,3E20.12)

		READ(10,'(A40)') ENDOFSECTION
		WRITE(9,'(A40)') ENDOFSECTION
		IF (ENDOFSECTION.EQ.'ENDOFSECTION') THEN
			WRITE(*,'(I6)') -1
		ELSE
			WRITE(*,*) 'Reading error, exit 2411'
			
			stop
		ENDIF
		
	!	Read elements
		READ(10,'(A30)') CONTRINFO
		WRITE(9,'(A30)') CONTRINFO
		NE=1
		DO WHILE(NE.LT.nFaces)
			READ(10,160) NE,NTYPE,NDP,(facePoints(NE,I),I=1,NDP)
			write(9,160) NE,NTYPE,NDP,(facePoints(NE,I),I=1,NDP)
			facePoints(NE,NDP+1)=facePoints(NE,1)
			if(NTYPE.eq.2) then     !Quadrilateral
				faceEdgesNum(NE)=4
			elseif(NTYPE.eq.3) then !triangle
				faceEdgesNum(NE)=3
			else
				write(*,*) 'Unsupported element types!'
				
				stop
			end if
		END DO
		160 Format(I8,1X,I2,1X,I2,1X,7I8:/(15X,7I8:))

		READ(10,'(A30)') ENDOFSECTION
		WRITE(9,'(A30)') ENDOFSECTION
		IF (ENDOFSECTION.EQ.'ENDOFSECTION') THEN
			WRITE(20,'(I6)') -1
		ELSE
			WRITE(*,*) 'Reading error, exit 2412'
		ENDIF

	!	Read element group control information record
		READ(10,'(A30)') CONTRINFO	
		WRITE(9,'(A30)') CONTRINFO	
		READ(10,*) GROUP,NGP,ELEMENT,NELGP,MATERIAL,MTYP,NSFLAGS,NFLAGS
		WRITE(9,*) GROUP,NGP,ELEMENT,NELGP,MATERIAL,MTYP,NSFLAGS,NFLAGS
		170 FORMAT(A8,I10,A11,I10,A10,I10,A8,I10)
		171 FORMAT(A8,1X, I10,1X,A11,1X,I10,1X,A10,1X,I10,1X,A8,1X,I10)
		180 FORMAT(A8,I10,A11,I10,A10,I10,A8,I10)

		READ(10,'(A40)') ELMAT
		WRITE(9,'(A40)') ELMAT
		READ(10,190) (ISOLVE(I),I=1,NFLAGS)
		WRITE(9,190) (ISOLVE(I),I=1,NFLAGS)
		190 FORMAT(10I8)
		READ(10,190) (NELT(I),I=1,NELGP)
		WRITE(9,190) (NELT(I),I=1,NELGP)

		READ(10,'(A30)') ENDOFSECTION
		WRITE(9,'(A30)') ENDOFSECTION

	!	Read boundary condition information records (total NBSETS)
	!   Only valid when IType=1 (element_side format)
		nBoundaryEdges=0
		do i=1,NBSETS
			READ(10,'(A30)') CONTRINFO	
			WRITE(9,'(A30)') CONTRINFO	
			READ(10,'(A32,8I10)') boundaryName, IType, NEntry, NValues, boundaryType
			WRITE(9,'(A32,8I10)') boundaryName, IType, NEntry, NValues, boundaryType
			if(IType.eq.0)then
				write(*,*) 'Only valid when the boundary is in element_side format!'
				
				stop
			end if
			
			do j=1,NEntry
				read(10,'(I10,I5,I5)') bFaceTemp, bFaceTypeTemp, bFaceDirTemp
				WRITE(9,'(I10,I5,I5)') bFaceTemp, bFaceTypeTemp, bFaceDirTemp
				boundaryEdgeMarkers(nBoundaryEdges+j)=boundaryType
				
				if(bFaceTypeTemp.eq.2) then !Quadrilateral
					pTemp1=facePoints(bFaceTemp,quadInfo(1,bFaceDirTemp))
					pTemp2=facePoints(bFaceTemp,quadInfo(2,bFaceDirTemp))
				elseif(bFaceTypeTemp.eq.3) then !Triangle
					pTemp1=facePoints(bFaceTemp,triInfo(1,bFaceDirTemp))
					pTemp2=facePoints(bFaceTemp,triInfo(2,bFaceDirTemp))
				else
					write(*,*) 'Unsupported element Types!'
					
					stop
				end if

				boundaryEdgePoints(nBoundaryEdges+j,1)=pTemp1
				boundaryEdgePoints(nBoundaryEdges+j,2)=pTemp2
			end do	

			nBoundaryEdges=nBoundaryEdges+NEntry

			READ(10,'(A30)') ENDOFSECTION
			WRITE(9,'(A30)') ENDOFSECTION
			IF (ENDOFSECTION.EQ.'ENDOFSECTION') THEN
				WRITE(*,'(I6)') -1
			ENDIF

		end do

		CLOSE(10)
		CLOSE(9)


	endif

	call cpu_time(time2)

	WRITE(*,*) 'CPU time = ', time2-time1


	end program DEMInterpolation

