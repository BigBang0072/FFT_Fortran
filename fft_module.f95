module FFT
contains

recursive subroutine FFT1D(total_length,total_depth,depth,start_index,omega,dump_credential,signal_array,transform_array)
    implicit none
    integer,intent(in) :: total_length,total_depth,depth,start_index
    real*4,dimension(0:total_length-1),intent(in) :: signal_array
    integer :: stride_length
    integer,dimension(1:total_depth,3),intent(inout) :: dump_credential
    integer :: start1,end1,start2,end2,i,j,counter !for iterating and combining step
    complex,dimension(0:total_length-1),intent(inout) :: transform_array
    complex :: temp_up,temp_down
    complex,intent(in) :: omega
    
    stride_length=2**depth    
    
    if(depth==total_depth+1)then
        !print *,"In last depth : ",dump_credential(depth-1,2) 
        !print *,"Start Index : ",start_index
        transform_array(dump_credential(depth-1,2))=cmplx(signal_array(start_index)*1,0.0) !Here 1 is FFT for 1 dimension.
    else
        !Calculating the first half of half-Fourier transform ,again will go recursive.
        !print *," ################# Start depth :",depth," #######################"
        !print *,"depth : ",depth
        dump_credential(depth,1)=1
        if(depth==1) then
            dump_credential(depth,2)=0
            start1=dump_credential(depth,2)
            dump_credential(depth,3)=(total_length/stride_length)-1
            end1=dump_credential(depth,3)
        else
            dump_credential(depth,2)=dump_credential(depth-1,2)
            start1=dump_credential(depth,2)
            dump_credential(depth,3)=dump_credential(depth-1,2)+(total_length/stride_length)-1
            end1=dump_credential(depth,3)
        end if
        !print *,"flag : ",dump_credential(depth,1),"starting : ",start1,&
                !"ending : ",end1
        call FFT1D(total_length,total_depth,depth+1,start_index,omega**2,&
                    dump_credential,signal_array,transform_array)
        
        !Calculating the secong half of the fourier transform.
        dump_credential(depth,1)=-1
        if(depth==1) then
            dump_credential(depth,2)=total_length/stride_length
            start2=dump_credential(depth,2)
            dump_credential(depth,3)=total_length-1
            end2=dump_credential(depth,3)
        else
            dump_credential(depth,2)=dump_credential(depth-1,2)+total_length/stride_length
            start2=dump_credential(depth,2)
            dump_credential(depth,3)=dump_credential(depth-1,3)
            end2=dump_credential(depth,3)
        end if
        !print *,"flag : ",dump_credential(depth,1),"starting : ",start2,&
                !"ending : ",end2
        call FFT1D(total_length,total_depth,depth+1,start_index+(stride_length/2),&
                    omega**2,dump_credential,signal_array,transform_array)
        
        !Combining the two half to backtrack our recursion.
        
        
        !do i=0,total_length-1
            !print *,transform_array(i)
        !end do
        
        counter=1
        !print *,"start1,end1,start2,end2 ",start1,end1,start2,end2
        !print *,"Omega**",omega
        
        do i=start1,end1
            j=end1+counter
            temp_up=transform_array(i)
            temp_down=transform_array(j)
            !Omega to be used is of previous level. eg.
            ! in last recombination ie FTn/2 and FTn/2 we use omega of FTn
            transform_array(i)=temp_up+(omega**(counter-1))*temp_down
            transform_array(j)=temp_up-(omega**(counter-1))*temp_down
            counter=counter+1
        end do
        
        !do i=0,total_length-1
            !print *,transform_array(i)
        !end do
        
        !print *,"###############  End Depth :",depth,"###########################"
    end if
    
end subroutine FFT1D

integer function log2(x)
    implicit none
    real,intent(in)::x
    
    log2=int(log(x)/log(2.0))
    
end function log2

end module FFT