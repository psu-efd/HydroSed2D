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


!   Read mesh from Gambit neutral format file
	subroutine readGambitMesh()
	USE COMMON_MODULE
	IMPLICIT NONE
	CHARACTER*80 NEUFILE,UNVFILE,CONTRINFO,GFILE,HED,PROG,REC5,&
			     ENDOFSECTION,MTHYEAR,GROUP,ELEMENT,MATERIAL,NSFLAGS,&
			     ELMAT,TODAY,MATNAM,boundaryName
	INTEGER NGRPS,NBSETS,NDFCD,NDFVL,ND,NE,i,j,&
			NTYPE,NDP,NODE(30),NGP,NELGP,MTYP,NFLAGS,ISOLVE(10),&
			NELT(1000000),REDATE(8),MATNUM,NODET,IType,NEntry,NValues,boundaryType
	integer bFaceTemp, bFaceTypeTemp, bFaceDirTemp,pTemp1,pTemp2
	real*8  boundarFaceTemp 
	REAL(8) REVL

	!information of edge for faces
	integer triInfo(2,3),quadInfo(2,4)
	data triInfo/1,2,2,3,3,1/,quadInfo/1,2,2,3,3,4,4,1/	

    !read data files from Gambit neutral format file of face
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
	READ(10,100) GFILE
	WRITE(*,100) GFILE
	READ(10,110) HED
	WRITE(*,110) HED
	READ(10,110) PROG
	WRITE(*,110) PROG
	READ(10,110) MTHYEAR
	READ(10,110) REC5
	READ(10,130) nNodes,nFaces,NGRPS,NBSETS,NDFCD,NDFVL
	READ(10,100) ENDOFSECTION
	100 FORMAT(A40)
	110 FORMAT(A80)
	120 FORMAT(A44,F5.2)
	130 FORMAT(6I10)
	140 FORMAT(5X,'nNodes',5X,'nFaces',5X,'NGRPS',4X,'NBSETS',5X,'NDFCD',5X,'NDFVL')
	
!	Read and write nodal coordinates
	READ(10,'(A30)') CONTRINFO
	WRITE(*,'(A30)') CONTRINFO
	ND=1
	DO WHILE(ND.LT.nNodes)
		READ(10,150) ND,pcoor(ND,1),pcoor(ND,2),pcoor(ND,3)
	END DO
	150 FORMAT(I10,1X,3E20.12)

	READ(10,'(A40)') ENDOFSECTION
	IF (ENDOFSECTION.EQ.'ENDOFSECTION') THEN
		WRITE(*,'(I6)') -1
	ELSE
		WRITE(*,*) 'Reading error, exit 2411'
		
		stop
	ENDIF
	
!	Read elements
	READ(10,'(A30)') CONTRINFO
	NE=1
	DO WHILE(NE.LT.nFaces)
		READ(10,160) NE,NTYPE,NDP,(facePoints(NE,I),I=1,NDP)
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
	IF (ENDOFSECTION.EQ.'ENDOFSECTION') THEN
		WRITE(20,'(I6)') -1
	ELSE
		WRITE(*,*) 'Reading error, exit 2412'
	ENDIF

!	Read element group control information record
	READ(10,'(A30)') CONTRINFO	
	READ(10,170) GROUP,NGP,ELEMENT,NELGP,MATERIAL,MTYP,NSFLAGS,NFLAGS
	170 FORMAT(A8,I10,A11,I10,A10,I10,A8,I10)
	180 FORMAT(A8,I10,A11,I10,A10,I10,A8,I10)

	READ(10,'(A40)') ELMAT
	READ(10,190) (ISOLVE(I),I=1,NFLAGS)
	190 FORMAT(10I8)
	READ(10,190) (NELT(I),I=1,NELGP)

	READ(10,'(A30)') ENDOFSECTION
	IF (ENDOFSECTION.EQ.'ENDOFSECTION') THEN
		WRITE(*,'(I6)') -1
	ENDIF

!	Read boundary condition information records (total NBSETS)
!   Only valid when IType=1 (element_side format)
	nBoundaryEdges=0
	do i=1,NBSETS
		READ(10,'(A30)') CONTRINFO	
		READ(10,'(A32,8I10)') boundaryName, IType, NEntry, NValues, boundaryType
		if(IType.eq.0)then
			write(*,*) 'Only valid when the boundary is in element_side format!'
			
			stop
		end if
		
		do j=1,NEntry
			read(10,'(I10,I5,I5)') bFaceTemp, bFaceTypeTemp, bFaceDirTemp
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
		IF (ENDOFSECTION.EQ.'ENDOFSECTION') THEN
			WRITE(*,'(I6)') -1
		ENDIF

	end do

	CLOSE(10)
	
	END SUBROUTINE
