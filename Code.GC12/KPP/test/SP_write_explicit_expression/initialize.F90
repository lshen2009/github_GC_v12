module initialize
implicit none
contains

  subroutine initialize_1D(x) 
  USE gckpp_Parameters
  implicit none
     REAL(kind=dp):: x(:)
     integer :: i
     CALL random_seed()
     do i=1,size(x)
       x(i)=rand()
     end do
  end subroutine initialize_1D

  subroutine num2str(num,strname)
    character(len=20) ::strname
    character(len=20)::string
    integer::num
    if (num<10) then
      string="(I1)"
    else if (num>=10 .and. num<100) then
      string="(I2)"
    else if (num>=100 .and. num<1000) then
      string="(I3)"
    else if (num>=1000 .and. num<10000) then
      string="(I4)"
    else if (num>=10000 .and. num<100000) then
      string="(I5)"
    endif
    write(strname,string) num
    strname=trim(strname)
  end subroutine num2str
  
end module initialize
