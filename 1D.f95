program fft_test
    use fft
    implicit none
    integer :: i,total_depth
    integer,parameter :: length=128*128*128
    integer,allocatable,dimension(:,:) :: dump_credential
    real*4,dimension(0:(length-1)) :: signal
    real*4,parameter :: pi=4*atan(1.0),tau=2*pi
    complex,parameter :: omega=cmplx(cos(tau/length),-sin(tau/length))
    complex,dimension(0:length-1) :: transform
    
    !real,external :: log2
    total_depth=log2(real(length))
    allocate(dump_credential(total_depth,3))
    !Generating the signal
    do i=0,length-1
        signal(i)=sin(i*tau/length)
        transform(i)=cmplx(0.0,0.0)
        print *,signal(i)
    end do
    
    !initializing Orientation Flag
    do i=1,total_depth
        dump_credential(i,:)=0
    end do
    
    !We will need seperate transform array.
    call FFT1D(length,total_depth,1,0,omega,dump_credential,signal,transform)
    
    do i=0,length-1
        print *,transform(i)
    end do
    
    deallocate(dump_credential)
   
end program fft_test


