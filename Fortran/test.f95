program test
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !! Perform unittests of the AQUA-Equation of state module
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
    implicit none

    integer,dimension(7) :: status_arr

    

    call test_LOAD_TABLE_PT(status_arr(1))
    call test_LOAD_TABLE_RhoT(status_arr(2))
    call test_LOAD_TABLE_RhoU(status_arr(3))

    print*,'------------------------------------------------------------'
    print*,'Test results:'
    print*,'------------------------------------------------------------'
    if (status_arr(1).eq.0)then
        print*,'    test_LOAD_TABLE_PT():          Passed'
    else
        print*,'    test_LOAD_TABLE_PT():          Failed, error:',status_arr(1)
    endif

    if (status_arr(2).eq.0)then
        print*,'    test_LOAD_TABLE_RhoT():        Passed'
    else
        print*,'    test_LOAD_TABLE_RhoT():        Failed, error:',status_arr(2)
    endif

    if (status_arr(3).eq.0)then
        print*,'    test_LOAD_TABLE_RhoU():        Passed'
    else
        print*,'    test_LOAD_TABLE_RhoU():        Failed, error:',status_arr(3)
    endif

    call test_BILINEAR_INTERPOLATION(status_arr(4))
    if (status_arr(4).eq.0)then
        print*,'    test_BILINEAR_INTERPOLATION(): Passed'
    else
        print*,'    test_BILINEAR_INTERPOLATION(): Failed, error:',status_arr(4)
    endif

    call test_INTERPOLATE_AQUA_PT(status_arr(5))
    if (status_arr(5).eq.0)then
        print*,'    test_INTERPOLATE_AQUA_PT():    Passed'
    else
        print*,'    test_INTERPOLATE_AQUA_PT():    Failed, error:',status_arr(5)
    endif

    call test_INTERPOLATE_AQUA_RhoT(status_arr(6))
    if (status_arr(6).eq.0)then
        print*,'    test_INTERPOLATE_AQUA_RhoT():  Passed'
    else
        print*,'    test_INTERPOLATE_AQUA_RhoT():  Failed, error:',status_arr(6)
    endif
print*,'------------------------------------------------------------'
if (any(status_arr.gt.0)) stop 1
end program test


subroutine test_LOAD_TABLE_PT(status)
    use aqua_eos,only:LOAD_TABLE_PT
    implicit none

    integer, intent(out) :: status
    CHARACTER(LEN=200) :: path

    status = 1

    path = '../Tables'
    CALL LOAD_TABLE_PT(path)
    
    ! if table could be loaded test was successfull
    status = 0
    return
end subroutine test_LOAD_TABLE_PT

subroutine test_LOAD_TABLE_RhoT(status)
    use aqua_eos,only:LOAD_TABLE_RhoT
    implicit none

    integer, intent(out) :: status
    CHARACTER(LEN=200) :: path

    status = 1

    path = '../Tables'
    CALL LOAD_TABLE_RhoT(path)
    
    ! if table could be loaded test was successfull
    status = 0
    return
end subroutine test_LOAD_TABLE_RhoT

subroutine test_LOAD_TABLE_RhoU(status)
    use aqua_eos,only:LOAD_TABLE_RhoU
    implicit none

    integer, intent(out) :: status
    CHARACTER(LEN=200) :: path

    status = 1

    path = '../Tables'
    CALL LOAD_TABLE_RhoU(path)
    
    ! if table could be loaded test was successfull
    status = 0
    return
end subroutine test_LOAD_TABLE_RhoU

subroutine test_BILINEAR_INTERPOLATION(status)
    use aqua_eos,only:BILINEAR_INTERPOLATION
    implicit none

    integer, intent(out) :: status
    !REAL(8),DIMENSION(size(zarr,1))     :: z
    real(8) :: x,x_arr(3)
    real(8) :: y,y_arr(3)
    real(8) :: z_arr2(2,3,3),z2(2),z_true
    status = 0

    x = 0.5_8
    y = 0.5_8
    x_arr = [0._8,1._8,2._8]
    y_arr = [0._8,1._8,2._8]
    z_arr2(:,:,1) = transpose(reshape([0.1_8,0.2_8,0.3_8,0.2_8,0.3_8,0.4_8],[3,2]))
    z_arr2(:,:,2) = transpose(reshape([0.2_8,0.3_8,0.4_8,0.2_8,0.3_8,0.4_8],[3,2]))
    z_arr2(:,:,3) = transpose(reshape([0.3_8,0.4_8,0.5_8,0.2_8,0.3_8,0.4_8],[3,2]))

    z2 = BILINEAR_INTERPOLATION(x,y,x_arr,y_arr,z_arr2)
    z_true = 0.2_8
    if (abs(z_true - z2(1)).gt. 1E-6_8 ) status = 1 
    z_true = 0.25_8
    if (abs(z_true - z2(2)).gt. 1E-6_8 ) status = 2 

    x = 1.0_8
    y = 0.25_8
    z2 = BILINEAR_INTERPOLATION(x,y,x_arr,y_arr,z_arr2)
    z_true = 0.225_8
    if (abs(z_true - z2(1)).gt. 1E-6_8 ) status = 3 
    z_true = 0.225_8
    if (abs(z_true - z2(2)).gt. 1E-6_8 ) status = 4 

    x = 1.5_8
    y = 1._8
    z2 = BILINEAR_INTERPOLATION(x,y,x_arr,y_arr,z_arr2)
    z_true = 0.35_8
    if (abs(z_true - z2(1)).gt. 1E-6_8 ) status = 3 
    z_true = 0.3_8
    if (abs(z_true - z2(2)).gt. 1E-6_8 ) status = 4 
    
    x = 2._8
    y = 1._8
    z2 = BILINEAR_INTERPOLATION(x,y,x_arr,y_arr,z_arr2)
    z_true = 0.4_8
    if (abs(z_true - z2(1)).gt. 1E-6_8 ) status = 5 
    z_true = 0.3_8
    if (abs(z_true - z2(2)).gt. 1E-6_8 ) status = 6 
    return
end subroutine test_BILINEAR_INTERPOLATION

subroutine test_INTERPOLATE_AQUA_PT(status)
    use aqua_eos,only: INTERPOLATE_AQUA_PT,safe_mode
    implicit none

    integer, intent(out) :: status
    real(8) :: p,t
    real(8),DIMENSION(8) :: EOS1,EOS2 !! [rho [kg/m^3], ad_grad [1], s [J/(kg K), u [J/kg], c [m/s]]
 
    status = 0

    t = 1000._8
    p = 1.E9_8

    safe_mode =.true.
    EOS1 = INTERPOLATE_AQUA_PT(P,T)
    safe_mode = .false.
    EOS2 = INTERPOLATE_AQUA_PT(P,T)

    if (any(EOS1.ne.EOS2)) status = 1
    
    t = 1.200113e+03_8
    p = 9.869633E+10_8
    safe_mode =.true.
    
    ! check against pure French & Redmer ice VII/X
    EOS2 = [3.090841e+03_8, 1.512685e-01_8,3.502381e+07_8*1e-4_8,-8.192848e+7_8+92.21378777E+6_8,&
     sqrt(1.494342e+8_8),18.01527E-3_8,0.0_8,0.0_8]
    EOS1 = INTERPOLATE_AQUA_PT(P,T)
    
    if (any(abs(EOS1-EOS2)/abs(EOS2+1E-20_8).gt.1E-3)) status =2
    
    return
end subroutine test_INTERPOLATE_AQUA_PT

subroutine test_INTERPOLATE_AQUA_RhoT(status)
    use aqua_eos,only: INTERPOLATE_AQUA_RhoT,safe_mode
    implicit none

    integer, intent(out) :: status
    real(8) :: rho,t
    real(8),DIMENSION(8) :: EOS1,EOS2 !! [rho [kg/m^3], ad_grad [1], s [J/(kg K), u [J/kg], c [m/s]]
 
    status = 0

    t = 1000._8
    rho = 911.75_8

    safe_mode =.true.
    EOS1 = INTERPOLATE_AQUA_RhoT(rho,T)
    safe_mode = .false.
    EOS2 = INTERPOLATE_AQUA_RhoT(rho,T)

    if (any(EOS1.ne.EOS2)) status = 1
    
    ! is differents since we use mazevet instead of French et al.
    t   = 8.02081497E+02_8
    rho = 2.09946287e+03_8
    safe_mode =.true.
    
    ! check against pure French & Redmer ice VII/X
    EOS2 = [2.20851640E+10_8, 9.92296970E-05*2.20851640E+10_8/(3.29585316E+03*rho),3.38718314e+07_8*1e-4_8,&
     -8.93316587e+7_8+92.21378777E+6_8,sqrt(4.30350277e+7_8),18.01527E-3_8,0.0_8,0.0_8]
    EOS1 = INTERPOLATE_AQUA_RhoT(rho,T)
    if (any(abs(EOS1-EOS2)/abs(EOS2+1E-20_8).gt.1E-3)) status = 2
    
    return
end subroutine test_INTERPOLATE_AQUA_RhoT

!! Leave that for creating new tests...
subroutine test_boilerplate(status)
    implicit none

    integer, intent(out) :: status

    status = 0
    return
end subroutine test_boilerplate