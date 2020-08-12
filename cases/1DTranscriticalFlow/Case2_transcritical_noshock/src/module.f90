MODULE  COMMON_MODULE 
		!plat form indicator: 1-windows(cvf6.0), 2-linux(gfortran)
		integer*4 platform
		parameter(platform = 1)

        integer*4 maxnodes_, maxedges_, maxfaces_, maxboundaryedges_,&
				  maxboundarypoints_,maxnodefaces_
        parameter(maxnodes_ = 15000, maxedges_ = 40000, maxfaces_ = 30000, &
		          maxboundaryedges_ = 40000,maxboundarypoints_=5000,maxnodefaces_=10)

		integer*4 ELEDGES       !maximum edges for each cell
		parameter(ELEDGES=3)

		INTEGER,PARAMETER::KKK=5
		!KKK:how many groups of sediment

        real*8 VSMALL, VLARGE
        parameter(VSMALL = 1e-32, VLARGE = 1e32)

		integer*4 meshType
		!mesh type: 1. gmsh  2. gambit neutral

		integer*4 curFace, curEdge, curPoint

!       Definition of data variables
!       1. points coordinats
        real*8 pcoor( maxnodes_, 3)        
!       points shared by triangle faces
        integer*4 pointFaces( maxnodes_, maxnodefaces_ )
!       how many faces are sharing this point
		integer*4 pointNFaces(maxnodes_)
		
!       2. edges (point IDs)
        integer*4 edgePoints( maxedges_, 2)
		integer*4 boundaryEdgePoints( maxboundaryedges_, 2)
		integer*4 boundaryEdgeGhostCells( maxfaces_ )
		integer*4 ghostCellsEdge( maxboundaryedges_ )
		integer*4 ghostCellsNeighbor (maxboundaryedges_ )
		real*8    edgeCenterCoor( maxedges_, 3)

!       3. faces (edge IDs, point IDs, only triangle for now)
		real*8    faceCenters( maxfaces_, 3)
        integer*4 faceEdges( maxfaces_, ELEDGES)
		integer*4 faceEdgesNum( maxfaces_ )
        integer*4 facePoints( maxfaces_, ELEDGES+1)
		integer*4 faceNeighbors( maxfaces_, ELEDGES)
		!faceNeighbors: neighbor cell number in three directions
		!faceNeighbors(i,j)=-1 if j direction is inlet
		!faceNeighbors(i,j)=-2 if j direction is outlet
		!faceNeighbors(i,j)=-3 if j direction is wall
				
		REAL(8),DIMENSION(maxfaces_,3)::faceLimiters
		!slope limiters for three conservative variables
		!see Anastasiou & Chan, 1997
		REAL(8),DIMENSION(maxfaces_)::etaLimiter
	
!		faces' edge normal vectors (pointing outside of the triangle)
		real*8 faceEdgeNormals( maxfaces_, ELEDGES, 3)

	    INTEGER,DIMENSION(maxfaces_, ELEDGES)::binfo
		!binfo:boundary info for each edge, 1-inlet, 2-outlet, 3-wall, 4-internal

		integer,dimension(maxfaces_)::pinfo
		!pinfo: position of the cell  1->inlet  2->outlet 3->walls 4->internal

!       4. topological data
!       edge's faces: owner and neighbour
!                     lower number face always be owner
!                     at least owner exists
        integer*4 edgeFaces( maxedges_, 2)

!       point location marker(where is the point located)
!       e.g.: 1    internal, 0   outside of domain
        integer*4 pointMarkers( maxnodes_ ) 

!       edge location marker(where is the edge located)
!       e.g.: 0    internal
!             1    wall
!             2    inlet
!             3    outlet
        integer*4 edgeMarkers( maxedges_ )
		integer*4 boundaryEdgeMarkers( maxboundaryedges_ )

!       edge length
		real*8 edgeLength( maxedges_ )

        integer*4 nNodes, nEdges, nFaces, nBoundaryEdges
		!nFaces: number of cells in the whole domain (wet plus dry)

	    integer*4 calc_flag( maxfaces_ ) 
!		calculation flag = 0 not in the calculation domain
!                        = 1 in the calculation domain

		REAL(8),DIMENSION(maxfaces_)::HiVi1,HiVi2,HiVi3,HiVi3_c,FdotN1,FdotN2,&
									  FdotN3,FdotN4_c,Q1,Q2,Q3,Q4_c
		!Q1:cell center xi value (free surface elevation)
		!Q2:cell center value of uh
		!Q3:cell center value of vh
		!FdotN1(maxfaces_)      = Result of vector operation F.n for continuity
		!FdotN2(maxfaces_)      = Result of vector operation F.n for x-momentum
		!FdotN3(maxfaces_)      = Result of vector operation F.n for y-momentum
		!HiVi1(maxfaces_)       = Source term vector component 1 - continuity
		!HiVi2(maxfaces_)       = Source term vector component 2 - x-momentum
		!HiVi3(maxfaces_)       = Source term vector component 3 - y-momentum
		!Q4_c:cell center value of sh
		!FdotN4_c: Result of vector operation F.n for sediment
		!Hivi4_c:  Source term vector component 3- sh

		REAL(8),DIMENSION(maxfaces_,2)::gradQ1,gradQ2,gradQ3,gradU,gradV,gradEta
		!gradQ1: cell center gradient of xi
		!gradQ2: cell center gradient of uh
		!gradQ3: cell center gradient of vh
		!gradU : cell center gradient of u
		!gradV : cell center gradient of v
		!gradEta : cell center gradient of eta

		REAL(8),DIMENSION(maxedges_,2)::edgeGradU,edgeGradV
		REAL(8),DIMENSION(maxedges_)::edgeQ1

		REAL(8),DIMENSION(maxfaces_)::Rem1,Rem2,Rem3,oldRem1,oldRem2,oldRem3
		!Rem1(maxfaces_)        \     this term = -F.n + Source terms
		!Rem2(maxfaces_)         >=   eq 4.41 in First Year report = d(Vq)/dt
		!Rem3(maxfaces_)        /     eq (22) in resubmitted paper = d(Vq)/dt
		!oldRem1(maxfaces_)     \
		!oldRem2(maxfaces_)      >=   Same as above but previous time step values
		!oldRem3(maxfaces_)     /     for Adams-Bashforth 2nd-order method (AB2)

		REAL(8),DIMENSION(maxfaces_)::Sox,Soy,vort,nb,coarse,UM,VN,eta
		!Sox(maxfaces_)         = Bottom slope in x-direction = -d(ZB)/dx
		!Soy(maxfaces_)         = Bottom slope in y-direction = -d(ZB)/dy
		!vort(maxfaces_)        = Depth-averaged Vorticity of computational cells
		!UM,VN:             cell center velocities
		!eta:               cell cneter water elevation (=h+zb)



		REAL(8),DIMENSION(maxfaces_)::DYY

		REAL(8),DIMENSION(maxfaces_,ELEDGES)::Qb1,Qb2,Qb3,Etab,lambda1,lambda2,lambda3
		!Qb1:Q1 at edge
		!Qb2:Q2 at edge
		!Qb3:Q3 at edge
		!Etab: eta at edge
		!lambda1(maxfaces_,3)   = Eigenvalue1 based on Roe-type averages
		!lambda2(maxfaces_,3)   = Eigenvalue1 based on Roe-type averages
		!lambda3(maxfaces_,3)   = Eigenvalue1 based on Roe-type averages

		!for ghost cells
		REAL(8),DIMENSION(maxboundaryedges_)::gQb1,gQb2,gQb3
		REAL(8),DIMENSION(maxboundaryedges_)::gFIofQb1,gFIofQb2,gFIofQb3
		REAL(8),DIMENSION(maxboundaryedges_)::gU,gV,gC,gAlpha
		REAL(8),DIMENSION(maxboundaryedges_)::gQ1,gQ2,gQ3,gZB,gEta

		REAL(8),DIMENSION(maxfaces_,ELEDGES)::FV1,FV2,FV3,FI1,FI2,FI3,FIofQb1,FIofQb2,FIofQb3
		REAL(8),DIMENSION(maxfaces_,ELEDGES)::u,v,alpha,Qb1av
		!u     = Roe-type Riemann average of u-velocity at cell interface   (Equation 16)
		!v     = Roe-type Riemann average of v-velocity at cell interface
		!alpha = Roe-type Riemann average of celerity at cell interface
		!Qb1av(maxfaces_,3)     = Roe-type Riemann average free-surface elevation at cell interface
		!FIofQb1(maxfaces_,3)   = Inviscid flux based on Riemann states at cell interface - continuity
		!FIofQb2(maxfaces_,3)   = Inviscid flux based on Riemann states at cell interface - x-momentum
		!FIofQb3(maxfaces_,3)   = Inviscid flux based on Riemann states at cell interface - y-momentum
		!FI1(maxfaces_,3)       = Inviscid flux based on Roe's flux funtion at cell interface - continuity 
		!FI2(maxfaces_,3)       = Inviscid flux based on Roe's flux funtion at cell interface - x-momentum
		!FI3(maxfaces_,3)       = Inviscid flux based on Roe's flux funtion at cell interface - y-momentum
		!FV1(maxfaces_,3)       = Viscous flux at cell interface based on cell centre values - continuity 
		!FV2(maxfaces_,3)       = Viscous flux at cell interface based on cell centre values - x-momentum
		!FV3(maxfaces_,3)       = Viscous flux at cell interface based on cell centre values - y-momentum

		REAL(8),DIMENSION(maxfaces_,ELEDGES,3,3)::R,L,modA,modLAMBDA
		!R(maxfaces_,3,3,3)     = Right Eigenvector matrix
		!L(maxfaces_,3,3,3)     = Left Eigenvector matrix
		!modLAMBDA(maxfaces_,3,3,3) = Matrix of absolute eigenvalues
		!modA(maxfaces_,3,3,3)  = Absolute flux Jabobian
		
		integer*4 runtype,inctr,outctr,frctl,sedimentctl,scalc,tscount,nDeltaT_output,&
		          nDeltaT_restart,sedInterval,sedTransRate 
		!inctr:  1-given h, 2-given u, 3-given q
		!outctr: 1-given h, 2-given u, 3-far field (zero gradient)
		!sedInterval: how many water time steps between sediment steps
		!sedTransRate: sediment transport rate formula 1-MPM, 2-Grass

		real*8 GrassConst(2) !constant for Grass sediment transport (A and m)

		real(8) inletH,inletUB,inletVB,inletQ,outletH,outletUB,outletVB,outletQ
		!inlet and outlet boundary conditions (not all are useful, depending on BC types)

		integer*4 inputcellno,hdoutputcellno
		REAL(8),DIMENSION(1000)::inputcellq,hdoutputcellq
		INTEGER,DIMENSION(1000)::inputcelli,hdoutputcelli
		!inputcelli:inlet cell numbers
		!hdoutputcelli:outlet cell numbers

		INTEGER::arrctl,terraindeal
		!terraindeal:地形处理,0网格高程采用单个值,1网格高程采用周边网格平均值
		!arrctl:选择恢复饱和系数的计算方法,0常规方法,1韦直林方法


		REAL(8)::vol,t,dt,tstop,g,a,b,d,visc,ts,td
		!a,b                = constants for integration using AB2
		!t                  = total time
		!dt                 = time step
		!visc               = kinematic eddy viscosuty coefficient
		!g                  = acceleration due to gravity 9.807m2/s
		!ts:hydrodynamics starting time
		!td:sediment starting time
		!tstop:stoping time
		!vol                = volume of total fluid in computational domain, based on cell centre variables
		
		REAL(8)::beta,tauwx,tauwy,Swx,Swy,Q2in
		!beta:limiter parameter, beta =1 -> minmod, beta=2 -> superbee
		!tauwx,tauwy        = wind stresses in x and y-directions, the steady state are set in subroutine 'wind_calc'
		!Swx,Swy            = Surface wind stress source terms  Swx=tauwx/rho, Swy=tauwy/rho
		!Q2in               = inlet VELOCITY, steady-state value set in subroutine 'uin_calc'

		REAL(8)::waterinput,wateroutput,watervolumn0

		REAL(8)::drydeep,mindeep,zaolv
		!drydeep:threshhold value for dry bed
		!mindeep:a value smaller than drydeep
		
		REAL(8)::watertemp1,watertemp2,watertemp3,sedtemp1,sedtemp2,sedtemp3
		REAL(8)::sedweight0,sedinput,sedoutput,termax,termin,hdoutq

		INTEGER::hp,fcp,fcpold,ltemp,div,qcp,qcells,mmax,lfno,arbcn
		INTEGER::ja,La,lfpn,nextlev,qncell,dcnter,done,fin,qqs,qqt,cklev
		
!       Sediment and other physical properties
!       rhos: sediment density   rhow: water density
!       diam: sediment diameter  porosity: sand porosity
!       thita_cri_b: critical Shields number for bed load
!       thita_cri_s: critical Shields number for suspended load
!       angleofRepose: angle of repose for sediment
!       g: gravity constant
!       C_slope: constant for slope effect
!       sigmac: Schmit number
!       Rs: submerged relative gravity of sediment
        real*8 rhos, rhow, diam, porosity, thita_cri_b, thita_cri_s, &
               angleofRepose, C_slope, sigmac, Rs

!       Shear stress vector at each node 
        real*8 nodeShearStress(maxnodes_, 2) !x,y direction components
!       Shear stress vector at each face centre 
        real*8 faceShearStress(maxfaces_, 2) !x,y direction components

!       face normal and steepest slope vector, face area projected on horizontal plane
        real*8 faceSteepestSlope(maxfaces_, 3)
        real*8 faceSlopeAngle(maxfaces_)
        real*8 face2DArea(maxfaces_)
		real*8 faceSubArea(maxfaces_, ELEDGES)
		!area of the triangle composed by face center and each three edges

!       Shields number and critial Shields number for each face after
!       slope effect correction
        real*8 faceShieldsNumber(maxfaces_)
        real*8 faceCriticalShieldsNumber(maxfaces_)

!       face center bed load sediment flux vector (qi)
        real*8 faceBedLoadFlux(maxfaces_, 2) !(now it is 2D vector)
!       edge center bed load sediment flus vector (qi) interpolated from face center
        real*8 edgeBedLoadFlux(maxedges_, 2)
!       divergence of bed load flux at face center (div(qi))
        real*8 divFaceBedLoadFlux(maxfaces_)
!       bed elevation change at face center
        real*8 faceElevationChange(maxfaces_)
!       bed elevation change at node
        real*8 nodeElevationChange(maxnodes_)

!		face Courant number
		real*8 faceCourantNumber(maxfaces_)
!		face nearest distance from face center to vertices
		real*8 faceMinR(maxfaces_)

!       max and min of Courant number in the domain
		real*8 maxCourantNumber,minCourantNumber

        real*8 qflat

		!variables for node values
		REAL(8),DIMENSION(maxnodes_)::nodeQ1,nodeQ2,nodeQ3,nodeZSurf
		REAL(8),DIMENSION(maxnodes_)::nodeU,nodeV,nodeSox,nodeSoy 
		REAL(8) meanH,meanU,meanV
        
        integer*4 ios, ierror

		character*25::meshfilename

		CHARACTER(len=100)::dirname,tempfilename,reloadfile
		integer*4 nStep
	    integer*4 ISTAT
	 
END MODULE  COMMON_MODULE 
