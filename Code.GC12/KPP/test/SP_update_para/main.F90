program main
  USE gckpp_LinearAlgebra
  IMPLICIT NONE

  INTEGER :: all_ind(234)
  LOGICAL :: species_ind(234)
  INTEGER, ALLOCATABLE :: new(:),x(:),y(:)
  INTEGER ::k,num
  DO k=1,234
   all_ind(k)=k
  ENDDO
  species_ind=(all_ind>230)
  new=pack(all_ind,species_ind)
  num=20
  ALLOCATE(x(num))

end program
