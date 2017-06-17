program test
    implicit none
    integer :: i,j
    integer,dimension(3,3) :: main
    integer,dimension(3) :: temp
    
    do i=1,3
        do j=1,3
            main(i,j)=i-j
        end do
        temp(i)=i
    end do
    
    do i=1,3
        do j=1,3
            print *,main(i,j)
        end do
        print *,"##"
        print *,temp(i)
        print *,"##"
    end do
    print *,"#####"
    temp(:)=main(:,1)
    do i=1,3
        print *,temp(i)
    end do
end program test