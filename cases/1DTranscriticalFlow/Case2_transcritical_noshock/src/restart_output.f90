!c***************************************************************
	subroutine restart_output
	USE COMMON_MODULE
	implicit none
 
	integer i,j,k,m
    CHARACTER*10 STRING

	if(mod(nStep,nDeltaT_restart).eq.0) then
		WRITE(UNIT=STRING, FMT='(I5)') int(nStep/nDeltaT_output)

	    if(platform.eq.1) then  !windows platform
		   tempfilename=trim(trim(trim(dirname) // '\results\re' // adjustl(STRING)) // '.dat')
	    else                    !linux platform 
	       tempfilename=trim(trim(trim(dirname) // '/results/re' // adjustl(STRING)) // '.dat')
	    endif

		open(unit=9,file=tempfilename,status='unknown')

 	    write(9,*) 'Reload file for SWEFVMSed'

		write(9,*) 'nNodes, nEdges, nFaces, nBoundaryEdges'
		write(9,*) nNodes, nEdges, nFaces, nBoundaryEdges

		write(9,*) 'curFace, curEdge, curPoint'
		write(9,*) curFace, curEdge, curPoint
		
		write(9,*) 'pcoor'
		write(9,*) ((pcoor(i,j),j=1,3),i=1,nNodes)

		write(9,*) 'pointFaces'
		write(9,*) ((pointFaces(i,j),j=1,maxnodefaces_),i=1,nNodes)

		write(9,*) 'pointNFaces'
		write(9,*) (pointNFaces(i),i=1,nNodes)		

		write(9,*) 'edgePoints'
		write(9,*) ((edgePoints(i,j),j=1,2),i=1,nEdges)

		write(9,*) 'boundaryEdgePoints'
		write(9,*) ((boundaryEdgePoints(i,j),j=1,2),i=1,maxboundaryedges_)

		write(9,*) 'boundaryEdgeGhostCells'
		write(9,*) (boundaryEdgeGhostCells(i),i=1,maxfaces_)		

		write(9,*) 'ghostCellsEdge'
		write(9,*) (ghostCellsEdge(i),i=1,maxboundaryedges_)		

		write(9,*) 'ghostCellsNeighbor'
		write(9,*) (ghostCellsNeighbor(i),i=1,maxboundaryedges_)		

		write(9,*) 'edgeCenterCoor'
		write(9,*) ((edgeCenterCoor(i,j),j=1,3),i=1,nEdges)

		write(9,*) 'faceCenters'
		write(9,*) ((faceCenters(i,j),j=1,3),i=1,nFaces)

		write(9,*) 'faceMinR'
		write(9,*) (faceMinR(i),i=1,nFaces)

		write(9,*) 'faceEdges'
		write(9,*) ((faceEdges(i,j),j=1,ELEDGES),i=1,nFaces)

		write(9,*) 'faceEdgesNum'
		write(9,*) (faceEdgesNum(i),i=1,nFaces)		

		write(9,*) 'facePoints'
		write(9,*) ((facePoints(i,j),j=1,ELEDGES+1),i=1,nFaces)

		write(9,*) 'faceNeighbors'
		write(9,*) ((faceNeighbors(i,j),j=1,ELEDGES),i=1,nFaces)

  		write(9,*) 'faceEdgeNormals'
		write(9,*) (((faceEdgeNormals(i,j,k),k=1,3),j=1,ELEDGES),i=1,nFaces)

  		write(9,*) 'binfo'
		write(9,*) ((binfo(i,j),j=1,ELEDGES),i=1,nFaces)

		write(9,*) 'pinfo'
		write(9,*) (pinfo(i),i=1,nFaces)		

  		write(9,*) 'edgeFaces'
		write(9,*) ((edgeFaces(i,j),j=1,2),i=1,nEdges)

		write(9,*) 'pointMarkers'
		write(9,*) (pointMarkers(i),i=1,nNodes)		

		write(9,*) 'edgeMarkers'
		write(9,*) (edgeMarkers(i),i=1,nEdges)		

		write(9,*) 'boundaryEdgeMarkers'
		write(9,*) (boundaryEdgeMarkers(i),i=1,maxboundaryedges_)		

		write(9,*) 'edgeLength'
		write(9,*) (edgeLength(i),i=1,nEdges)		

		write(9,*) 'edgeQ1'
		write(9,*) (edgeQ1(i),i=1,nEdges)		

		write(9,*) 'calc_flag'
		write(9,*) (calc_flag(i),i=1,nFaces)		
				
		write(9,*) 'Q1'
		write(9,*) (Q1(i),i=1,nFaces)

		write(9,*) 'Q2'
		write(9,*) (Q2(i),i=1,nFaces)

		write(9,*) 'Q3'
		write(9,*) (Q3(i),i=1,nFaces)

		write(9,*) 'Q4_c'
		write(9,*) (Q4_c(i),i=1,nFaces)

		write(9,*) 'Rem1'
		write(9,*) (Rem1(i),i=1,nFaces)
		write(9,*) 'Rem2'
		write(9,*) (Rem2(i),i=1,nFaces)
		write(9,*) 'Rem3'
		write(9,*) (Rem3(i),i=1,nFaces)
		write(9,*) 'oldRem1'
		write(9,*) (oldRem1(i),i=1,nFaces)
		write(9,*) 'oldRem2'
		write(9,*) (oldRem2(i),i=1,nFaces)
		write(9,*) 'oldRem3'
		write(9,*) (oldRem3(i),i=1,nFaces)

		write(9,*) 'Sox'
		write(9,*) (Sox(i),i=1,nFaces)
		write(9,*) 'Soy'
		write(9,*) (Soy(i),i=1,nFaces)
		write(9,*) 'vort'
		write(9,*) (vort(i),i=1,nFaces)
		write(9,*) 'nb'
		write(9,*) (nb(i),i=1,nFaces)
		write(9,*) 'coarse'
		write(9,*) (coarse(i),i=1,nFaces)
		write(9,*) 'UM'
		write(9,*) (UM(i),i=1,nFaces)
		write(9,*) 'VN'
		write(9,*) (VN(i),i=1,nFaces)
		write(9,*) 'eta'
		write(9,*) (eta(i),i=1,nFaces)

		write(9,*) 'DYY'
		write(9,*) (DYY(i),i=1,nFaces)
		
		write(9,*) 'gQ1'
		write(9,*) (gQ1(i),i=1,maxboundaryedges_)
		write(9,*) 'gQ2'
		write(9,*) (gQ2(i),i=1,maxboundaryedges_)
		write(9,*) 'gQ3'
		write(9,*) (gQ3(i),i=1,maxboundaryedges_)
		write(9,*) 'gZB'
		write(9,*) (gZB(i),i=1,maxboundaryedges_)
		write(9,*) 'gEta'
		write(9,*) (gEta(i),i=1,maxboundaryedges_)

		write(9,*) 'inputcellno,hdoutputcellno'
		write(9,*) inputcellno,hdoutputcellno

		write(9,*) 'inputcellq,hdoutputcellq'
		write(9,*) (inputcellq(i),i=1,1000)
		write(9,*) (hdoutputcellq(i),i=1,1000)

		write(9,*) 'inputcelli,hdoutputcelli'
		write(9,*) (inputcelli(i),i=1,1000)
		write(9,*) (hdoutputcelli(i),i=1,1000)

		write(9,*) 'waterinput,wateroutput,watervolumn0'
		write(9,*) waterinput,wateroutput,watervolumn0
		
		write(9,*) 'watertemp1,watertemp2,watertemp3,sedtemp1,sedtemp2,sedtemp3'
		write(9,*) watertemp1,watertemp2,watertemp3,sedtemp1,sedtemp2,sedtemp3
		write(9,*) 'sedweight0,sedinput,sedoutput,termax,termin,hdoutq'
		write(9,*) sedweight0,sedinput,sedoutput,termax,termin,hdoutq
	
  		write(9,*) 'faceSteepestSlope'
		write(9,*) ((faceSteepestSlope(i,j),j=1,3),i=1,nFaces)

		write(9,*) 'faceSlopeAngle'
		write(9,*) (faceSlopeAngle(i),i=1,nFaces)

		write(9,*) 'face2DArea'
		write(9,*) (face2DArea(i),i=1,nFaces)

  		write(9,*) 'faceSubArea'
		write(9,*) ((faceSubArea(i,j),j=1,ELEDGES),i=1,nFaces)

		write(9,*) 'meshfilename'
		write(9,*) meshfilename

        write(9,*) 'dirname'
		write(9,*) dirname
		
        write(9,*) 'tempfilename'
		write(9,*) tempfilename

        write(9,*) 'reloadfile'
		write(9,*) reloadfile
		
		close(9)
		
	end if

	end subroutine
