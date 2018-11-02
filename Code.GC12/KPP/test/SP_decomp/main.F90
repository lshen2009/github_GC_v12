program main
  IMPLICIT NONE

  INTEGER :: all_ind(234)
  LOGICAL :: ind(234)
  INTEGER, ALLOCATABLE :: new(:)
  INTEGER ::k
  DO k=1,234
   all_ind(k)=k
  ENDDO
  ind=(all_ind>230)
  new=pack(all_ind,ind)
  print *, new
  print *, sum(array=all_ind,mask=.not.ind)
end program
