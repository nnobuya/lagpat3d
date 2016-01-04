module mod_cnst

  implicit none

  private
  public:: npt, ndim, pi, rm_sol, npt_x1, npt_x2, npt_x3

  !..part. num. and dimension
  integer, parameter:: ndim = 3
  !integer, parameter:: npt_x1 = 40, npt_x2 = 40, npt_x3 = 20
  integer, parameter:: npt_x1 = 20, npt_x2 = 20, npt_x3 = 20
  integer, parameter:: npt = npt_x1 *npt_x2 *npt_x3

  !integer, parameter:: nin = 1, nou = 2  !! for hokan grid


  !..const.
  real(8), parameter:: pi     = 3.141592653589793d0
  real(8), parameter:: rm_sol = 1.9891d33

end module mod_cnst
