! =================================================================
! Generic comments about the structure of the program and short HOWTO's
! =================================================================

program main
  use constants
  use params
  ! use utils, only:
  ! OTHER program-level MODULES
  implicit none

  ! Arrays allocated based on the inputs
  real(wp),allocatable :: array(:,:)   !> [Dim1 by Dim2] Unit?
  integer  :: i !> Dummy integer
  integer(li)  :: j !> Long integer
  real(wp) :: var !> Variable description here

  call read_parameters()

  ! INITIALIZE THE CALCULATION
  write(6,*) "Information  about running the program:           "!, some_information
  write(6,*) ! Use blank lines

  ! RUN THE PROGRAM
  write(6,*)
  write(6,*) "RUN: "  
  ! call run(parameters)

end program main