program fft_test
    use fft
    implicit none
    integer :: i,j,k
    integer,parameter :: dim1=4,dim2=4,dim3=4
    real*4,parameter :: pi=4*atan(1.0),tau=2*pi
    complex,dimension(0:dim1-1,0:dim2-1,0:dim3-1) :: signal
    complex,dimension(0:dim1-1,0:dim2-1,0:dim3-1) :: transform
    
    
    !Generating Signal
    do i=0,dim1-1
        do j=0,dim2-1
            do k=0,dim3-1
                signal(i,j,k)=cmplx(sin(i*tau/dim1+j*tau/dim2+k*tau/dim3),0.0)
                print *,signal(i,j,k)
                transform(i,j,k)=cmplx(0.0,0.0)
            end do
        end do
    end do
    
    !Calling Transform
    call FFT3D(dim1,dim2,dim3,signal,transform)
    
    print *,"#####################################"
    !Printing Transform
    do i=0,dim1-1
        do j=0,dim2-1
            do k=0,dim3-1
                print *,transform(i,j,k)
            end do
            print *,"   #####"
        end do
        print *,"##############"
    end do
    
end program fft_test