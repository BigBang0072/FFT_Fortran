program fft_test
    use fft
    implicit none
    integer :: i,j
    integer,parameter :: dim1=8,dim2=8
    complex,dimension(0:dim1-1,0:dim2-1) :: signal
    real*4,parameter :: pi=4*atan(1.0),tau=2*pi
    complex,parameter :: omega1=cmplx(cos(tau/dim1),-sin(tau/dim1))
    complex,parameter :: omega2=cmplx(cos(tau/dim2),-sin(tau/dim2))
    complex,dimension(0:dim1-1,0:dim2-1) :: transform
    
    !Generating the Signal
    do i=0,dim1-1
        do j=0,dim2-1
            signal(i,j)=cmplx(sin(i*(tau/dim1)+j*(tau/dim2)),0.0)
            transform(i,j)=cmplx(0.0,0.0)
        end do
    end do
    
    call FFT2D(dim1,dim2,signal,transform)
    
    !Printing the values
    do i=0,dim1-1
        do j=0,dim2-1
            print *,transform(i,j)
        end do
        print *,"#####################"
    end do
    
end program fft_test