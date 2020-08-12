!c***************************************************************
	subroutine results_output
	USE COMMON_MODULE,ONLY: Q1,Q2,Q3,dirname,tempfilename,nStep,&
							faceCenters,nDeltaT_output,eta,UM,VN,nodeU,&
							nodeV,facePoints,nNodes,nFaces,pcoor,pointMarkers,&
							meanH,meanU,meanV,nodeZSurf,faceEdgesNum,ELEDGES,faceLimiters,&
							nodeElevationChange,Sox,Soy,nodeSox,nodeSoy,platform,nodeQ1
	implicit none
 
	integer i,j
    CHARACTER*10 STRING

	if((mod(nStep,nDeltaT_output).eq.0).or.(nStep.eq.0)) then
		WRITE(UNIT=STRING, FMT='(I5)') int(nStep/nDeltaT_output)

		!write cell center variables
!		tempfilename=trim(trim(trim(dirname) // '\results\phi' // adjustl(STRING)) // '.dat')
!		open(unit=9,file=tempfilename,status='unknown')

!		write(9,*) 'VARIABLES = "X", "Y", "Z", "phi1","phi2","phi3"'

!		do i=1,nFaces
!			write(9,*) faceCenters(i,1),faceCenters(i,2),faceCenters(i,2),&
!					   faceLimiters(i,1),faceLimiters(i,2),faceLimiters(i,3)
!		end do
!		close(9)

		!interpolate cell center data to nodes
		call cellCenterToNodes(Sox,nodeSox)
		call cellCenterToNodes(Soy,nodeSoy)
		call cellCenterToNodes(eta,nodeZSurf)

		!calculate mean H, U, and V
		call meanHUV

		!write water surface elevation file(for node values)
	    if(platform.eq.1) then  !windows platform
		    tempfilename=trim(trim(trim(dirname) // '\results\h' // adjustl(STRING)) // 'tec.dat')
	    else                    !linux platform 
	        tempfilename=trim(trim(trim(dirname) // '/results/h' // adjustl(STRING)) // 'tec.dat')
	    endif

		open(unit=9,file=tempfilename,status='unknown')

 	    write(9,*) 'TITLE = SWE Solver Water Surface Elevation'
		write(9,*) 'VARIABLES = "X", "Y", "Z", "Zsurf","H", "U", "V", "Umag","deltaZ","Sox","Soy"'

		if(ELEDGES.eq.3) then
			write(9, *) 'ZONE, DATAPACKING=POINT, N=', nNodes, ',E=', &
				         nFaces,' ,ZONETYPE=FETRIANGLE'
		elseif(ELEDGES.eq.4) then
			write(9, *) 'ZONE, DATAPACKING=POINT, N=', nNodes, ',E=', &
				         nFaces, ' ,ZONETYPE=FEQUADRILATERAL'
		end if

		do i=1,nNodes
			if(pointMarkers(i).eq.1) then   !internal point
				write(9, *) pcoor(i,1), pcoor(i,2), pcoor(i,3), nodeZSurf(i),nodeQ1(i),&
					        nodeU(i), nodeV(i), dsqrt(nodeU(i)**2+nodeV(i)**2), nodeElevationChange(i),&
							nodeSox(i),nodeSoy(i)
			else                            !external point
				write(9, *) pcoor(i,1), pcoor(i,2), pcoor(i,3), nodeZSurf(i),meanH,&
					        meanU, meanV, dsqrt(meanU**2+meanV**2), ' 0.0 ', '0 ','0'
			end if

		enddo

	    do i=1,nFaces
			write(9, *) (facePoints(i,j),j=1,faceEdgesNum(i))
		enddo

		close(9)
		
	end if

	end subroutine
