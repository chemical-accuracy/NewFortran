module params
  use constants
  implicit none

integer,parameter :: MAXPARAM=1000 !> Some maximum parameter
logical :: FLAGS(MAXPARAM) !> can be used to change behavior of program
real(wp) :: parameter !> parameter that will be read
contains

  subroutine read_parameters()
    read(5,*)                        ! skip header line
    read(5,*)  parameter             ! read in some parameter

    ! ! Get the job ID
    ! call getenv('SLURM_ARRAY_TASK_ID',JobID)
    ! JobID = adjustl(JobID)
    write(6,*) "INPUT:"
    write(6,*)
    write(6,*) "Parameter        :           ", parameter
  end subroutine read_parameters

end module params
