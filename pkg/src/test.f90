subroutine update_alpha(Y,n)
  implicit none
  integer :: Y(*)
  integer :: n
  print *,Y(1),n
  Y(1)=Y(1)+1
  print *, Y(1),n
end subroutine update_alpha
