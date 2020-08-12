MODULE  COMMON_MODULE 
		!plat form indicator: 1-windows(cvf6.0), 2-linux(gfortran)
		integer*4 platform
		parameter(platform = 2)

        integer*4 maxnodes_, maxedges_, maxfaces_, maxboundaryedges_,&
				  maxboundarypoints_,maxnodefaces_
        parameter(maxnodes_ = 211580, maxedges_ = 610000, maxfaces_ = 210000, &
		          maxboundaryedges_ = 10000,maxboundarypoints_=5000,maxnodefaces_=10)

		integer*4 ELEDGES       !maximum edges for each cell
		parameter(ELEDGES=4)

        real*8 VSMALL, VLARGE
        parameter(VSMALL = 1e-32, VLARGE = 1e32)

		integer*4 meshType
		!mesh type: 1. gmsh  2. gambit neutral

		integer*4 curFace, curEdge, curPoint

        integer*4 nNodes, nEdges, nFaces, nBoundaryEdges
		!nFaces: number of cells in the whole domain (wet plus dry)

		integer*4 nElements


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
!       e.g.: 0   outside of domain
!             1   inlet boundary
!             2   outlet boundary
!             3   wall boundary (has the top priority among three boundaries)
!             4    internal 
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

       
        integer*4 ios, ierror

		character*25::meshfilename

		CHARACTER(len=100)::dirname,tempfilename,reloadfile
		integer*4 nStep
	    integer*4 ISTAT
	 
END MODULE  COMMON_MODULE 
