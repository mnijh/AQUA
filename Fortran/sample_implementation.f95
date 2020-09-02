PROGRAM sample_implementation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! Sample implementation of the aqua_eos module.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! Copyright (c) 2020, Jonas Haldemann
    !
    ! Permission is hereby granted, free of charge, to any person obtaining a copy of this software    
    ! and associated documentation files (the "Software"), to deal in the Software without restriction,
    ! including without limitation the rights to use, copy, modify, merge, publish, distribute,        
    ! sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is    
    ! furnished to do so, subject to the following conditions:
    !
    ! The above copyright notice and this permission notice shall be included in all copies or         
    ! substantial portions of the Software.
    ! 
    ! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
    ! NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND           
    ! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,     
    ! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   
    ! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.'
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    USE AQUA_EOS,only:DIM_OUTPUT,LOAD_TABLE_PT,INTERPOLATE_AQUA_PT
    USE AQUA_EOS,only:LOAD_TABLE_RhoT,INTERPOLATE_AQUA_RHOT
    USE AQUA_EOS,only:LOAD_TABLE_RhoU,INTERPOLATE_AQUA_RHOU
    IMPLICIT NONE
    INTEGER  :: ii
    REAL(8)  :: T,P,lP,rho,U
    REAL(8)  :: EOS(DIM_OUTPUT)
    CHARACTER(LEN=200) :: path

    !! Load the AQUA tables
    path = '../Tables'
    CALL LOAD_TABLE_PT(path)
    CALL LOAD_TABLE_RhoT(path)
    CALL LOAD_TABLE_RhoU(path)

    T = 1400._8 ! [K]
    lP = 1._8  ! [Pa]
    
    do ii = 0,110
        P = 10._8**(lP + ii/10._8) 
	!! Evaluate the tabulated data
        EOS = INTERPOLATE_AQUA_PT(P,T)
        write(*,*)'*************************************************'
        write(*,*)'Pressure [Pa]:            ',P
        write(*,*)'Temperartue [K]:          ',T
        write(*,*)'Density [kg/m^3]:         ',EOS(1)
        write(*,*)'Adiabatic_grad  :         ',EOS(2)
        write(*,*)'Entropy [J/(kg K)]:       ',EOS(3)
        write(*,*)'Internal Energy [J/(kg)]: ',EOS(4)
        write(*,*)'Bulk speed of sound [m/s]:',EOS(5)
    enddo

END PROGRAM sample_implementation
