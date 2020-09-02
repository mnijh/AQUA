MODULE AQUA_EOS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! Fortran 95 module to evaluate the AQUA-Equation of state (EoS) Tables
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
    !
    ! If used for a scientific publication please cite:
    !     @article{Haldemann_2020a,
    !       doi = {},
    !       url = {},
    !       year  = {2020},
    !       month = {},
    !       volume = {},
    !       number = {},
    !       pages = {},
    !       author = {Haldemann, J., Alibert, Y., Mordasini, Ch. and Benz, W. },
    !       title = {AQUA: A Collection of $H_2O$ Equations of State for Planetary Models},
    !       journal = {Astronomy & Astrophysics}
    !     }
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    IMPLICIT NONE

    INTEGER,PARAMETER  :: WP = 8                                     !! Precision of real variables
    INTEGER,PARAMETER  :: un = 42                                    !! Input unit
    INTEGER,PARAMETER  :: DIM_OUTPUT = 8                             !! DIMENSION of output arrays
    CHARACTER(LEN=20)  :: table_name_pt   = 'aqua_eos_pt_v1_0.dat'   !! Name of PT_table
    CHARACTER(LEN=22)  :: table_name_rhot = 'aqua_eos_rhot_v1_0.dat' !! Name of RhoT_table
    CHARACTER(LEN=22)  :: table_name_rhou = 'aqua_eos_rhou_v1_0.dat' !! Name of RhoU_table

    REAL(WP),DIMENSION(:),ALLOCATABLE :: P_arr         !! Pressure range of PT-AQUA-table
    REAL(WP),DIMENSION(:),ALLOCATABLE :: T1_arr        !! Temperature range of PT-AQUA-table
    REAL(WP),DIMENSION(:),ALLOCATABLE :: T2_arr        !! Temperature range of RHOT-AQUA-table
    REAL(WP),DIMENSION(:),ALLOCATABLE :: Rho1_arr      !! Temperature range of RHOT-AQUA-table
    REAL(WP),DIMENSION(:),ALLOCATABLE :: Rho2_arr      !! Temperature range of RHOU-AQUA-table
    REAL(WP),DIMENSION(:),ALLOCATABLE :: U_arr         !! Temperature range of RHOU-AQUA-table
    REAL(WP),DIMENSION(:,:,:),ALLOCATABLE :: EOS_PT    !! AQUA-EoS data
    REAL(WP),DIMENSION(:,:,:),ALLOCATABLE :: EOS_RhoT  !! AQUA-EoS data
    REAL(WP),DIMENSION(:,:,:),ALLOCATABLE :: EOS_RhoU  !! AQUA-EoS data

    INTEGER :: nP   = 1093
    INTEGER :: nT   = 301
    INTEGER :: nRho = 1501
    INTEGER :: nU   = 461

    LOGICAL :: safe_mode            = .true.           !! Check if input is inside table boundaries
    LOGICAL :: notify_outside_table = .true.           !! Notify if outside table boundaries (only checked in safe_mode)
    LOGICAL :: stop_outside_table   = .true.           !! Stop if outside of table (only checked in safe_mode) 
                                                       !! otherwise evaluate at last P,T inside tabluated range. 
    
    CONTAINS

    SUBROUTINE LOAD_TABLE_PT(path)
        !! Load EoS table located in path. Check if file exists otherwise stop with error.
        IMPLICIT NONE
        
        INTEGER,PARAMETER :: skip_first_lines = 19 !! Skip first n lines of table

        CHARACTER(LEN=*),INTENT(IN) :: path        !! Path to directory with aqua table
        LOGICAL  :: iexist
        REAL(WP) :: ptmp,ttmp
        INTEGER  :: ii,jj,phase

        !! Display logo and license statement. Do not erease or comment this part.
        CALL LOGO
        CALL LICENSE

        !! Check if file exists
        inquire(file=trim(path)//'/'//trim(table_name_pt),exist=iexist)
        if (iexist) then
            write(*,*)'[AQUA-EOS] -- Loading Table: ',trim(path)//'/'//trim(table_name_pt)
            open(unit=un,file=trim(path)//'/'//trim(table_name_pt),action='read')
        else
            write(*,*)'ERROR: [AQUA-EOS] -- File: ',trim(path)//'/'//trim(table_name_pt),' does not exist!'
            STOP
        endif

        !! 
        ALLOCATE(P_arr(nP))
        ALLOCATE(T1_arr(nT))
        ALLOCATE(EOS_PT(DIM_OUTPUT,nT,nP))

        do ii=1,skip_first_lines
            READ(un,*)    
        enddo
        
        !! Read line by line 
        do ii=1,nP
            do jj=1,nT
                READ(un,*) ptmp,ttmp,EOS_PT(1:DIM_OUTPUT,jj,ii),phase
                P_arr(ii)  = log10(ptmp)
                T1_arr(jj) = log10(ttmp)
            enddo
        enddo

        write(*,*)'[AQUA-EOS] -- Finished Loading Table'

    END SUBROUTINE LOAD_TABLE_PT 
    
    SUBROUTINE LOAD_TABLE_RhoT(path)
        !! Load EoS table located in path. Check if file exists otherwise stop with error.
        IMPLICIT NONE
        
        INTEGER,PARAMETER :: skip_first_lines = 21 !! Skip first n lines of table

        CHARACTER(LEN=*),INTENT(IN) :: path        !! Path to directory with aqua table
        LOGICAL  :: iexist
        REAL(WP) :: rhotmp,ttmp
        INTEGER  :: ii,jj,phase

        !! Display logo and license statement. Do not erease or comment this part.
        CALL LOGO
        CALL LICENSE

        !! Check if file exists
        inquire(file=trim(path)//'/'//trim(table_name_rhot),exist=iexist)
        if (iexist) then
            write(*,*)'[AQUA-EOS] -- Loading Table: ',trim(path)//'/'//trim(table_name_rhot)
            open(unit=un,file=trim(path)//'/'//trim(table_name_rhot),action='read')
        else
            write(*,*)'ERROR: [AQUA-EOS] -- File: ',trim(path)//'/'//trim(table_name_rhot),' does not exist!'
            STOP
        endif

        !! 
        ALLOCATE(Rho1_arr(nRho))
        ALLOCATE(T2_arr(nT))
        ALLOCATE(EOS_RhoT(DIM_OUTPUT,nT,nRho))

        do ii=1,skip_first_lines
            READ(un,*)    
        enddo
        
        !! Read line by line 
        do ii=1,nRho
            do jj=1,nT
                READ(un,*) rhotmp,ttmp,EOS_RhoT(1:DIM_OUTPUT,jj,ii),phase
                Rho1_arr(ii) = log10(rhotmp)
                T2_arr(jj)   = log10(ttmp)
            enddo
        enddo

        write(*,*)'[AQUA-EOS] -- Finished Loading Table'

    END SUBROUTINE LOAD_TABLE_RhoT 

    SUBROUTINE LOAD_TABLE_RhoU(path)
        !! Load EoS table located in path. Check if file exists otherwise stop with error.
        IMPLICIT NONE
        
        INTEGER,PARAMETER :: skip_first_lines = 21 !! Skip first n lines of table

        CHARACTER(LEN=*),INTENT(IN) :: path        !! Path to directory with aqua table
        LOGICAL  :: iexist
        REAL(WP) :: rhotmp,utmp
        INTEGER  :: ii,jj,phase

        !! Display logo and license statement. Do not erease or comment this part.
        CALL LOGO
        CALL LICENSE

        !! Check if file exists
        inquire(file=trim(path)//'/'//trim(table_name_rhou),exist=iexist)
        if (iexist) then
            write(*,*)'[AQUA-EOS] -- Loading Table: ',trim(path)//'/'//trim(table_name_rhou)
            open(unit=un,file=trim(path)//'/'//trim(table_name_rhou),action='read')
        else
            write(*,*)'ERROR: [AQUA-EOS] -- File: ',trim(path)//'/'//trim(table_name_rhou),' does not exist!'
            STOP
        endif

        !! 
        ALLOCATE(Rho2_arr(nRho))
        ALLOCATE(U_arr(nU))
        ALLOCATE(EOS_RhoU(DIM_OUTPUT,nU,nRho))

        do ii=1,skip_first_lines
            READ(un,*)    
        enddo
        
        !! Read line by line 
        do ii=1,nRho
            do jj=1,nU
                READ(un,*) rhotmp,utmp,EOS_RhoU(1:DIM_OUTPUT,jj,ii),phase
                Rho2_arr(ii) = log10(rhotmp)
                U_arr(jj) = log10(utmp)
            enddo
        enddo

        write(*,*)'[AQUA-EOS] -- Finished Loading Table'

    END SUBROUTINE LOAD_TABLE_RhoU 

    FUNCTION INTERPOLATE_AQUA_PT(P,T) result(EOS)
        !! Interpolate AQUA-EoS table given pressure P [Pa] and temperature T [K].
        !! Returns state variables i.e. [rho [kg/m^3], ad_grad [1], s [J/(kg K), u [J/kg], c [m/s]]
        IMPLICIT NONE

        REAL(WP),INTENT(IN) :: P     !! Pressure [Pa]
        REAL(WP),INTENT(IN) :: T     !! Temperature [K]
        REAL(WP),DIMENSION(DIM_OUTPUT) :: EOS !! [rho [kg/m^3], ad_grad [1], s [J/(kg K), u [J/kg], c [m/s]]

        REAL(WP) :: lP,lT

        if (safe_mode) then
            if ((log10(P).ge.P_arr(1)).and.(log10(P).le.P_arr(nP))) then
                lP = log10(P)
            else
                if (notify_outside_table) then
                    write(*,*) 'WARNING: [AQUA-EOS (INTERPOLATE_AQUA_PT)] -- Pressure outside tabulated range.'
                    write(*,*) 'WARNING: [AQUA-EOS (INTERPOLATE_AQUA_PT)] -- Pressure:',P
                endif
                if (stop_outside_table) STOP
                if (log10(P).lt.P_arr(1)) then
                    lP = P_arr(1)
                else
                    lP = P_arr(nP)
                endif
            endif

            if ((log10(T).ge.T1_arr(1)).and.(log10(T).le.T1_arr(nT))) then
                lT = log10(T)
            else
                if (notify_outside_table) then
                    write(*,*) 'WARNING: [AQUA-EOS (INTERPOLATE_AQUA_PT)] -- Temperature outside tabulated range.'
                    write(*,*) 'WARNING: [AQUA-EOS (INTERPOLATE_AQUA_PT)] -- Temperature:',T
                endif
                if (stop_outside_table) STOP
                if (log10(T).lt.T1_arr(1)) then
                    lT = T1_arr(1)
                else
                    lT = T1_arr(nT)
                endif
            endif
        else
            lP = log10(P)
            lT = log10(T)
        endif

        EOS = BILINEAR_INTERPOLATION(lP,lT,P_arr,T1_arr,EOS_PT)

    END FUNCTION INTERPOLATE_AQUA_PT

    FUNCTION INTERPOLATE_AQUA_RhoT(Rho,T) result(EOS)
        !! Interpolate AQUA-EoS table given density Rho [kg/m^3] and temperature T [K].
        !! Returns state variables i.e. [pressure [Pa], ad_grad [1], s [J/(kg K), u [J/kg], c [m/s]]
        IMPLICIT NONE

        REAL(WP),INTENT(IN) :: Rho   !! Density [kg/m^3]
        REAL(WP),INTENT(IN) :: T     !! Temperature [K]
        REAL(WP),DIMENSION(DIM_OUTPUT) :: EOS !! [p [Pa], ad_grad [1], s [J/(kg K), u [J/kg], c [m/s]]

        REAL(WP) :: lRho,lT

        if (safe_mode) then
            if ((log10(Rho).ge.Rho1_arr(1)).and.(log10(Rho).le.Rho1_arr(nRho))) then
                lRho = log10(Rho)
            else
                if (notify_outside_table) then
                    write(*,*) 'WARNING: [AQUA-EOS (INTERPOLATE_AQUA_RhoT)] -- Density outside tabulated range.'
                    write(*,*) 'WARNING: [AQUA-EOS (INTERPOLATE_AQUA_RhoT)] -- Density:',Rho
                endif
                if (stop_outside_table) STOP
                if (log10(Rho).lt.Rho1_arr(1)) then
                    lRho = Rho1_arr(1)
                else
                    lRho = Rho1_arr(nRho)
                endif
            endif

            if ((log10(T).ge.T2_arr(1)).and.(log10(T).le.T2_arr(nT))) then
                lT = log10(T)
            else
                if (notify_outside_table) then
                    write(*,*) 'WARNING: [AQUA-EOS (INTERPOLATE_AQUA_RhoT)] -- Temperature outside tabulated range.'
                    write(*,*) 'WARNING: [AQUA-EOS (INTERPOLATE_AQUA_RhoT)] -- Temperature:',T
                endif
                if (stop_outside_table) STOP
                if (log10(T).lt.T2_arr(1)) then
                    lT = T2_arr(1)
                else
                    lT = T2_arr(nT)
                endif
            endif
        else
            lRho = log10(Rho)
            lT = log10(T)
        endif

        EOS = BILINEAR_INTERPOLATION(lRho,lT,Rho1_arr,T2_arr,EOS_RhoT)

    END FUNCTION INTERPOLATE_AQUA_RhoT

    FUNCTION INTERPOLATE_AQUA_RhoU(Rho,U) result(EOS)
        !! Interpolate AQUA-EoS table given density Rho [kg/m^3] and internal energy U [J/kg].
        !! Returns state variables i.e. [pressure [Pa], temperature [K], ad_grad [1], s [J/(kg K), c [m/s]]
        IMPLICIT NONE

        REAL(WP),INTENT(IN) :: Rho   !! Density [kg/m^3]
        REAL(WP),INTENT(IN) :: U     !! Internal Energy [J/kg]
        REAL(WP),DIMENSION(DIM_OUTPUT) :: EOS !! [p [Pa], t [K], ad_grad [1], s [J/(kg K), c [m/s]]

        REAL(WP) :: lRho,lU

        if (safe_mode) then
            if ((log10(Rho).ge.Rho2_arr(1)).and.(log10(Rho).le.Rho2_arr(nRho))) then
                lRho = log10(Rho)
            else
                if (notify_outside_table) then
                    write(*,*) 'WARNING: [AQUA-EOS (get_AQUA_EOS)] -- Density outside tabulated range.'
                    write(*,*) 'WARNING: [AQUA-EOS (get_AQUA_EOS)] -- Density:',Rho
                endif
                if (stop_outside_table) STOP
                if (log10(Rho).lt.Rho2_arr(1)) then
                    lRho = Rho2_arr(1)
                else
                    lRho = Rho2_arr(nRho)
                endif
            endif

            if ((log10(U).ge.U_arr(1)).and.(log10(U).le.U_arr(nU))) then
                lU = log10(U)
            else
                if (notify_outside_table) then
                    write(*,*) 'WARNING: [AQUA-EOS (get_AQUA_EOS)] -- Internal Energy outside tabulated range.'
                    write(*,*) 'WARNING: [AQUA-EOS (get_AQUA_EOS)] -- Internal Energy:',U
                endif
                if (stop_outside_table) STOP
                if (log10(U).lt.U_arr(1)) then
                    lU = U_arr(1)
                else
                    lU = U_arr(nU)
                endif
            endif
        else
            lRho = log10(Rho)
            lU = log10(U)
        endif

        EOS = BILINEAR_INTERPOLATION(lRho,lU,Rho2_arr,U_arr,EOS_RhoU)

    END FUNCTION INTERPOLATE_AQUA_RhoU

    FUNCTION BILINEAR_INTERPOLATION(x,y,xarr,yarr,zarr)
        !! Bilienear interpolaion of vector z(x,y).
        IMPLICIT NONE
        REAL(WP),INTENT(IN)                  :: x,y
        REAL(WP),INTENT(IN),DIMENSION(:)     :: xarr,yarr
        REAL(WP),INTENT(IN),DIMENSION(:,:,:) :: zarr          !! zarr = [nz,ny,nx]
        REAL(WP),DIMENSION(size(zarr,1))     :: BILINEAR_INTERPOLATION

        INTEGER   :: nx,ny,ix0,iy0,ix1,iy1,ii
        REAL(WP)  :: dx,dy
        
        ! get array length
        nx = size(xarr)
        ny = size(yarr)

        ! get array spacing
        dx = xarr(2)-xarr(1)
        dy = yarr(2)-yarr(1)

        ! get neighbors within array
        ix0 = max(1,min(1 + FLOOR((x-xarr(1))/dx),nx))
        if (mod(x-xarr(1),dx).eq.0._WP) then
            ix1 = ix0
        else
            ix1 = min(ix0 + 1,nx)
        endif

        ! get neighbors within array
        iy0 = max(1,min(1 + FLOOR((y-yarr(1))/dy),ny))
        if (mod(y-yarr(1),dy).eq.0._WP) then
            iy1 = iy0
        else
            iy1 = min(iy0 + 1,ny)
        endif
        
        ! check if x-array is not properly spaced.
        if (x.gt.xarr(ix1)) then
            do ii = 1,nx-ix1
                if (x.le.xarr(ix1+ii)) then
                    ix0=ix1+ii-1
                    ix1=ix0+1
                    EXIT
                endif
            enddo
        else if (x.lt.xarr(ix0)) then
            do ii = 1,ix0
                if (x.ge.xarr(ix0-ii)) then
                    ix0=ix0-ii
                    ix1=ix0+1
                    EXIT
                endif
            enddo
        endif

        ! check if y-array is not properly spaced.
        if (y.gt.yarr(iy1)) then
            do ii = 1,ny-iy1
                if (y.le.yarr(iy1+ii)) then
                    iy0=iy1+ii-1
                    iy1=iy0+1
                    EXIT
                endif
            enddo
        else if (y.lt.yarr(iy0)) then
            do ii = 1,iy0
                if (y.ge.yarr(iy0-ii)) then
                    iy0=iy0-ii
                    iy1=iy0+1
                    EXIT
                endif
            enddo
        endif

        ! check if input point is in one of the input arrays
        if (ix0.eq.ix1) then
            if (iy0.eq.iy1) then
                BILINEAR_INTERPOLATION = zarr(:,iy0,ix0)
            else
                BILINEAR_INTERPOLATION = interpol_1d(y,yarr(iy0),yarr(iy1),zarr(:,iy0,ix0),zarr(:,iy1,ix0))
            endif
        else if (iy0.eq.iy1) then
            BILINEAR_INTERPOLATION = interpol_1d(x,xarr(ix0),xarr(ix1),zarr(:,iy0,ix0),zarr(:,iy0,iy1))
        else
            BILINEAR_INTERPOLATION = interpol_2d(x,y,xarr(ix0),xarr(ix1),yarr(iy0),yarr(iy1),&
                zarr(:,iy0,ix0),zarr(:,iy0,ix1),zarr(:,iy1,ix0),zarr(:,iy1,ix1))
        endif

        CONTAINS

            PURE FUNCTION interpol_1d(x,x0,x1,y0,y1)
                !! Linear interpolation
                IMPLICIT NONE

                REAL(WP),INTENT(IN)              :: x,x0,x1
                REAL(WP),DIMENSION(:),INTENT(IN) :: y0,y1
                REAL(WP),DIMENSION(size(y0))     :: interpol_1d 

                interpol_1d = y0 + ( x - x0 ) * (y1 - y0) / (x1 - x0)
                
                return
            end FUNCTION interpol_1d

            PURE FUNCTION interpol_2d(x,y,x0,x1,y0,y1,z00,z10,z01,z11)
                !! Bilinear interpolation
                IMPLICIT NONE

                REAL(WP),INTENT(IN)              :: x,y,x0,x1,y0,y1
                REAL(WP),DIMENSION(:),INTENT(IN) :: z00,z10,z01,z11
                REAL(WP),DIMENSION(size(z00))    :: interpol_2d

                interpol_2d = ( z00*(x1-x)*(y1-y) + z10*(x-x0)*(y1-y)  &
                  + z01*(x1-x)*(y-y0) + z11*(x-x0)*(y-y0) ) / ( (x1-x0) * (y1-y0) )
                
                return
            end FUNCTION interpol_2d
    END FUNCTION BILINEAR_INTERPOLATION

    SUBROUTINE LICENSE()
        !! Write licences statement        
        write(*,*)'Copyright (c) 2020, Jonas Haldemann'
        write(*,*)''
        write(*,*)'Permission is hereby granted, free of charge, to any person obtaining a copy  '
        write(*,*)'of this software and associated documentation files (the "Software"), to deal '
        write(*,*)'in the Software without restriction, including without limitation the rights  '
        write(*,*)'to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     '
        write(*,*)'copies of the Software, and to permit persons to whom the Software is         '
        write(*,*)'furnished to do so, subject to the following conditions:'
        write(*,*)''
        write(*,*)'The above copyright notice and this permission notice shall be included in all'
        write(*,*)'copies or substantial portions of the Software.'
        write(*,*)''
        write(*,*)'THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    '
        write(*,*)'IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      '
        write(*,*)'FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   '
        write(*,*)'AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        '
        write(*,*)'LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, '
        write(*,*)'OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS        '
        write(*,*)'IN THE SOFTWARE.'
        write(*,*)''

    END SUBROUTINE LICENSE

    SUBROUTINE LOGO
        write(*,*)"                                                                             "
        write(*,*)"              .o.         .oooooo.      ooooo     ooo       .o.              "
        write(*,*)"             .888.       d8P'  `Y8b     `888'     `8'      .888.             "
        write(*,*)'            .8"888.     888      888     888       8      .8"888.            '
        write(*,*)"           .8' `888.    888      888     888       8     .8' `888.           "
        write(*,*)"          .88ooo8888.   888      888     888       8    .88ooo8888.          "
        write(*,*)"         .8'     `888.  `88b    d88b     `88.    .8'   .8'     `888.         "
        write(*,*)"        o88o     o8888o  `Y8bood8P'Ybd'    `YbodP'    o88o     o8888o        "
        write(*,*)"-----------------------------------------------------------------------------"
        write(*,*)'        A Collection of H2O Equations of States for Planetary Models         '
        write(*,*)"-----------------------------------------------------------------------------"
        write(*,*)''
    END SUBROUTINE LOGO

END MODULE AQUA_EOS
