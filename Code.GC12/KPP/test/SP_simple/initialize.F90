module initialize
implicit none
contains

  subroutine initialize_1D(x) 
  implicit none
     real:: x(:)
     integer :: i
     CALL random_seed()
     do i=1,size(x)
       x(i)=rand()
     end do
  end subroutine initialize_1D

end module initialize
