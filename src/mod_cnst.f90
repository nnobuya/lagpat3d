module mod_cnst

  implicit none

  private
  public:: npt, ndim, pi, rm_sol, npt_x1, npt_x2, npt_x3, set_cnst

  !..part. num. and dimension
  integer:: ndim, npt, npt_x1, npt_x2, npt_x3

  !..const.
  real(8), parameter:: pi     = 3.141592653589793d0
  real(8), parameter:: rm_sol = 1.9891d33

contains

  subroutine set_cnst

    implicit none


    !..part. num. and dimension
    ndim = 3
    npt_x1 = 20
    npt_x2 = 20
    npt_x3 = 20

    npt = npt_x1 *npt_x2 *npt_x3

    return

  end subroutine set_cnst

end module mod_cnst
