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

!   Initialization
    subroutine init
	USE COMMON_MODULE
	implicit none
    
	integer ii,jj,k,mm,o
	REAL*8 BEDPKK(5),EMPKK(5)
	character*20 temp,string
	character*14 tem
	character*15 terrainfile
	character*15 qtfile,qhfile,stfile

	integer GETCWD

	ISTAT = GETCWD(dirname)

!c      *********  opening 'ctr.dat'  ********************
	if(platform.eq.1) then  !windows platform
	     tempfilename=trim(trim(dirname) // '\input\ctr.dat')
	else                    !linux platform 
	     tempfilename=trim(trim(dirname) // '/input/ctr.dat')
	endif


	open(8, file = tempfilename, status = 'old', iostat = ios)
    
	read(8,*)string
	read(8,*)runtype
	read(8,*)string
	read(8,*)reloadfile
	read(8,*)string
	read(8,*)dt,sedInterval,ts,td,tstop
	read(8,*)string
	read(8,*)g,visc,beta,zaolv
	read(8,*)string
	read(8,*)drydeep,mindeep
	read(8,*)string
	read(8,*)inctr
	read(8,*)string
	read(8,*)inletH,inletUB,inletVB,inletQ	
	read(8,*)string
	read(8,*)outctr
	read(8,*)string
	read(8,*)outletH,outletUB,outletVB,outletQ
	read(8,*)string
	read(8,*)nDeltaT_output
	read(8,*)string
	read(8,*)nDeltaT_restart
	read(8,*)string
	read(8,*)frctl
	read(8,*)string
	read(8,*)sedimentctl
	read(8,*)string
	read(8,*)terraindeal,termax,termin
	read(8,*)string
	read(8,*)terrainfile              !terrainfile is a local variable
	read(8,*)string
	read(8,*)meshType
	read(8,*)string
	read(8,*)meshfilename
	read(8,*)string
	read(8,*)qtfile
	read(8,*)string
	read(8,*)qhfile
	read(8,*)string
	read(8,*)stfile
	close(unit=8)

!   read sediment properties
	call readSedimentProperties

	tscount=0

	if(runtype==0)then
		write(*,*) 'Run from start'

		nStep=0
			
		modLAMBDA=0D0
		R=0D0
		L=0D0
		oldRem1=0D0
		oldRem2=0D0
		oldRem3=0D0
		Qb1=0D0
		Qb2=0D0
		Qb3=0D0

		gQb1=0D0
		gQb2=0D0
		gQb3=0D0

		
		waterinput=0
		wateroutput=0
		calc_flag=0

		Swx=0.0
		Swy=0.0
		
		!read mesh file
		if(meshType.eq.1) then   !gmsh file
			call readGMSHMesh
		elseif(meshType.eq.2) then !gambit neutral file
			call readGambitMesh
		else
			write(*,*) 'Unsupported mesh format!'
			
			stop
		endif

		!build mesh information
		call buildMeshData

		call IC_setup
		t=0.0D0
	else
		write(*,*) 'Run from restart'
		call read_reload
		call IC_setup
		t=t+dt
	end if

	return
	end


!c  This subroutine sets the coefficients used in the integration
	subroutine int_coeffs
	USE COMMON_MODULE,ONLY: a,b,d,t,dt

	implicit none

	if(t.lt.2.0*dt)then
		a=1.0D0
		b=0.0D0
		d=0.0D0
	else
		a=3.0D0/2.0D0
		b=-1.0D0/2.0D0
	end if

	end subroutine

!c  This subroutine set up some initial conditions including:
!c  Q1(i), Q2(i), Q3(i)
!c  Sox(i), Soy(i)
!c  Vi(i), dl(i)
	subroutine IC_setup
	USE COMMON_MODULE,ONLY: nFaces,Q1,Q2,Q3,&
		Sox,Soy,waterinput,&
		wateroutput,watervolumn0,sedweight0,sedinput,sedoutput,&
		nb,frctl,sedimentctl,&
		binfo,Qb1,Qb2,Qb3,drydeep,&
		mindeep,zaolv,inctr,outctr,&
		pinfo,t,dt,ts,td,&
		face2DArea,faceCenters,nBoundaryEdges,&
		gQ1,gQ2,gQ3,gZB,ghostCellsNeighbor,inletH,inletQ,UM,VN,eta,&
		nodeQ1,nodeQ2,nodeQ3,nodeU,nodeV,nodeZSurf,gEta

	implicit none
	integer:: i,j,k,counttemp,rposx,rposy
	integer dirx(8),diry(8)
	real(8):: ttemp,initq
	integer NoUse

	watervolumn0=0
	sedweight0=0
	ttemp=0.0

!	call get_input_qtotal(ttemp,initq)
!	call from_q_get_h(initq,initz)
!   initz=0.5

	do i=1,nFaces	      
		nb(i)=zaolv

		Q1(i)=inletH-faceCenters(i,3)
	    Q2(i)=inletQ							 !uh on center
	    Q3(i)=0.0D0									 !vh on center

		!Hard wire some initial conditions
		!dam break for this case
!		if(faceCenters(i,1).le.2.39) then
!			Q1(i)=0.25+0.33-faceCenters(i,3)
!		end if

		if(Q1(i)<=drydeep)then	
			UM(i)=0
			VN(i)=0
			eta(i)=faceCenters(i,3)
		else
			UM(i)=Q2(i)/Q1(i)
			VN(i)=Q3(i)/Q1(i)
			eta(i)=faceCenters(i,3)+Q1(i)
		end if
	end do

	!for ghost cells, variables as gXXX
	do i=1,nBoundaryEdges
		if(frctl.eq.1)then   !complicated bottom
			gZB(i)=faceCenters(ghostCellsNeighbor(i),3)
		else                 !uniform bottom
			gZB(i)=0D0
		end if

		gEta(i)=eta(ghostCellsNeighbor(i))
		gQ1(i)=Q1(ghostCellsNeighbor(i))
		gQ2(i)=0D0
		gQ3(i)=0D0

	enddo

    !find all inlet cells, only for east or west inlet now
	call find_inputcells()

	do i=1,nFaces
		watervolumn0=watervolumn0+Q1(i)*face2DArea(i)
    end do       

  	!interpolate cell center Q, U and V to nodes
	call cellCenterToNodes(Q1,nodeQ1)
	call cellCenterToNodes(Q2,nodeQ2)
	call cellCenterToNodes(Q3,nodeQ3)
	call cellCenterToNodes(UM,nodeU)
	call cellCenterToNodes(VN,nodeV)
	call cellCenterToNodes(eta,nodeZsurf)


	if(sedimentctl.eq.1)then
	end if
	
	return
	
	end

!c*******************************************************************************
	subroutine IC_setup_reload
	USE COMMON_MODULE
	implicit none
	integer:: i,j,k


	return
	end


!c
!c********************************************************************
!c     find all inlet cells (cells with a boundary of inlet   
	subroutine find_inputcells()
	USE COMMON_MODULE, only:inputcellno,inputcelli,faceEdges,edgeMarkers,nFaces

	implicit none
	integer i

	inputcellno=0
	do i=1,nFaces
		if(edgeMarkers(faceEdges(i,1))==1.or. &
		   edgeMarkers(faceEdges(i,2))==1.or. &
		   edgeMarkers(faceEdges(i,3))==1) then
			inputcellno=inputcellno+1
			inputcelli(inputcellno)=i
		end if
	end do
	return
	end subroutine

!c
!c********************************************************************
!c  find all the outlet cells (cells with a boundary of outlet
	subroutine find_hd_outputcells()
	USE COMMON_MODULE,ONLY: hdoutputcellno,hdoutputcelli,faceEdges,edgeMarkers,nFaces
	implicit none
	integer i

	hdoutputcellno=0
	do i=1,nFaces
		if(edgeMarkers(faceEdges(i,1))==2.or. &
		   edgeMarkers(faceEdges(i,1))==2.or. &
		   edgeMarkers(faceEdges(i,1))==2) then
			hdoutputcellno=hdoutputcellno+1
			hdoutputcelli(hdoutputcellno)=i
		end if
	end do

	return
	end subroutine

!c***************************************************************
	subroutine read_reload
	USE COMMON_MODULE
	implicit none
 
	integer i,j,k,m
    CHARACTER*10 STRING,temp

	 if(platform.eq.1) then  !windows platform
	    tempfilename=trim(trim(trim(dirname) // '\results\' // reloadfile))
	 else                    !linux platform 
	    tempfilename=trim(trim(trim(dirname) // '/results/' // reloadfile))
	 endif

	open(unit=9,file=tempfilename,status='unknown', iostat = ios)
    if ( ios /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Could not open the reload file'
        
	    stop
    else
 	    read(9,*) temp

		read(9,*) temp
		read(9,*) nNodes, nEdges, nFaces, nBoundaryEdges

		read(9,*) temp
		read(9,*) curFace, curEdge, curPoint
		
		read(9,*) temp
		read(9,*) ((pcoor(i,j),j=1,3),i=1,nNodes)

		read(9,*) temp
		read(9,*) ((pointFaces(i,j),j=1,maxnodefaces_),i=1,nNodes)

		read(9,*) temp
		read(9,*) (pointNFaces(i),i=1,nNodes)		

		read(9,*) temp
		read(9,*) ((edgePoints(i,j),j=1,2),i=1,nEdges)

		read(9,*) temp
		read(9,*) ((boundaryEdgePoints(i,j),j=1,2),i=1,maxboundaryedges_)

		read(9,*) temp
		read(9,*) (boundaryEdgeGhostCells(i),i=1,maxfaces_)		

		read(9,*) temp
		read(9,*) (ghostCellsEdge(i),i=1,maxboundaryedges_)		

		read(9,*) temp
		read(9,*) (ghostCellsNeighbor(i),i=1,maxboundaryedges_)		

		read(9,*) temp
		read(9,*) ((edgeCenterCoor(i,j),j=1,3),i=1,nEdges)

		read(9,*) temp
		read(9,*) ((faceCenters(i,j),j=1,3),i=1,nFaces)

		read(9,*) temp
		read(9,*) (faceMinR(i),i=1,nFaces)

		read(9,*) temp
		read(9,*) ((faceEdges(i,j),j=1,ELEDGES),i=1,nFaces)

		read(9,*) temp
		read(9,*) (faceEdgesNum(i),i=1,nFaces)		

		read(9,*) temp
		read(9,*) ((facePoints(i,j),j=1,ELEDGES+1),i=1,nFaces)

		read(9,*) temp
		read(9,*) ((faceNeighbors(i,j),j=1,ELEDGES),i=1,nFaces)

  		read(9,*) temp
		read(9,*) (((faceEdgeNormals(i,j,k),k=1,3),j=1,ELEDGES),i=1,nFaces)

  		read(9,*) temp
		read(9,*) ((binfo(i,j),j=1,ELEDGES),i=1,nFaces)

		read(9,*) temp
		read(9,*) (pinfo(i),i=1,nFaces)		

  		read(9,*) temp
		read(9,*) ((edgeFaces(i,j),j=1,2),i=1,nEdges)

		read(9,*) temp
		read(9,*) (pointMarkers(i),i=1,nNodes)		

		read(9,*) temp
		read(9,*) (edgeMarkers(i),i=1,nEdges)		

		read(9,*) temp
		read(9,*) (boundaryEdgeMarkers(i),i=1,maxboundaryedges_)		

		read(9,*) temp
		read(9,*) (edgeLength(i),i=1,nEdges)	
		
		read(9,*) temp
		read(9,*) (edgeQ1(i),i=1,nEdges)		
	
		read(9,*) temp
		read(9,*) (calc_flag(i),i=1,nFaces)		
				
		read(9,*) temp
		read(9,*) (Q1(i),i=1,nFaces)

		read(9,*) temp
		read(9,*) (Q2(i),i=1,nFaces)

		read(9,*) temp
		read(9,*) (Q3(i),i=1,nFaces)

		read(9,*) temp
		read(9,*) (Q4_c(i),i=1,nFaces)

		read(9,*) temp
		read(9,*) (Rem1(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (Rem2(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (Rem3(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (oldRem1(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (oldRem2(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (oldRem3(i),i=1,nFaces)

		read(9,*) temp
		read(9,*) (Sox(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (Soy(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (vort(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (nb(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (coarse(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (UM(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (VN(i),i=1,nFaces)
		read(9,*) temp
		read(9,*) (eta(i),i=1,nFaces)

		read(9,*) temp
		read(9,*) (DYY(i),i=1,nFaces)
		
		read(9,*) temp
		read(9,*) (gQ1(i),i=1,maxboundaryedges_)
		read(9,*) temp
		read(9,*) (gQ2(i),i=1,maxboundaryedges_)
		read(9,*) temp
		read(9,*) (gQ3(i),i=1,maxboundaryedges_)
		read(9,*) temp
		read(9,*) (gZB(i),i=1,maxboundaryedges_)
		read(9,*) temp
		read(9,*) (gEta(i),i=1,maxboundaryedges_)

		read(9,*) temp
		read(9,*) inputcellno,hdoutputcellno

		read(9,*) temp
		read(9,*) (inputcellq(i),i=1,1000)
		read(9,*) (hdoutputcellq(i),i=1,1000)

		read(9,*) temp
		read(9,*) (inputcelli(i),i=1,1000)
		read(9,*) (hdoutputcelli(i),i=1,1000)
		
		read(9,*) temp
		read(9,*) waterinput,wateroutput,watervolumn0
		
		read(9,*) temp
		read(9,*) watertemp1,watertemp2,watertemp3,sedtemp1,sedtemp2,sedtemp3
		read(9,*) temp
		read(9,*) sedweight0,sedinput,sedoutput,termax,termin,hdoutq

  		read(9,*) temp
		read(9,*) ((faceSteepestSlope(i,j),j=1,3),i=1,nFaces)

		read(9,*) temp
		read(9,*) (faceSlopeAngle(i),i=1,nFaces)

		read(9,*) temp
		read(9,*) (face2DArea(i),i=1,nFaces)

  		read(9,*) temp
		read(9,*) ((faceSubArea(i,j),j=1,ELEDGES),i=1,nFaces)
        
		read(9,*) temp
		read(9,*) meshfilename

        read(9,*) temp
		read(9,*) dirname
		
        read(9,*) temp
		read(9,*) tempfilename

        read(9,*) temp
		read(9,*) reloadfile
		
		close(9)
		
	end if

	end subroutine
