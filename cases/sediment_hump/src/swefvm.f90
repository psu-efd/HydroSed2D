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

!   cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   c Finite Volume Method for the 2D Shallow Water Equation
!   c with Sediment Transport on Unstructured Mesh
!   c (the sediment transport part is solved seperately with 
!   c the flow.)
!   c by Xiaofeng Liu, Hydrosystems Lab, UIUC, 10/26/2008
!   C liu19@illinois.edu

    program main
	USE COMMON_MODULE
    implicit none
	real*8 timestart,timestop,timetotal

	write(*,*) '-------------------------------------------------------------------'
    write(*,*) '|   Shallow    | HydroSed2D: Shallow Water Equation with Sediment |'
	write(*,*) '|   Water      |                                                  |'
	write(*,*) '|   Equation   |    by     : Xiaofeng Liu, Hydrosystems Lab, UIUC |'
	write(*,*) '|   with       |                                                  |'
	write(*,*) '|   Sediment   |    version: 1.0  (12/5/2008)                     |' 
	write(*,*) '-------------------------------------------------------------------'
    write(*,*) 'Web  : http://vtchl.uiuc.edu/~liu19     '
	write(*,*) 'Email: liu19@illinois.edu'
	write(*,*) ' '
	
	call cpu_time(timestart)

	!initialization
	call init

	call results_output

	!calculation: main time step loop
	do while (t <= tstop) 
		write(*,*) 't = ', t, ' s out of ', tstop, 's' 
		nStep=nStep+1
		call swe
		
		if((sedimentctl.eq.1).and.(scalc.eq.1)) then
			call sediment
		end if
		
		call results_output
		call restart_output

	    t=t+dt

		tscount=tscount+1
		if(tscount==sedInterval)then
			scalc=1	
			tscount=0
		else
			scalc=0
		end if

    end do

	call cpu_time(timestop)

	timetotal = timestop - timestart

	write(*,*) 'Total CPU time = ', timetotal, ' seconds.'

	end
