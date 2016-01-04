module mod_set

  use mod_cnst, only: npt, ndim

  implicit none

  private
  public:: k_zoku, i_test, last_lp, int_t, int_x, nin, nou, &
       & nout_lpt, n_anim, &
       & x1_in, x1_out, x2_in, x2_out, x3_in, x3_out, &
       & set_param, nx1, nx2, nx3, &
       & d_fld, t_fld, ye_fld, en_fld, x_fld, dx_fld, dma, v_fld, v0_fld, &
       & set_data

  integer:: k_zoku, i_test, last_lp, int_t, int_x, nin, nou
  integer:: nout_lpt, n_anim
  real(8):: r_in, r_out, bound_in, bound_out
  real(8):: x1_in, x1_out, x2_in, x2_out, x3_in, x3_out

  !..main
  integer:: nx1, nx2, nx3

  !..Grid & field data (hydro results)

  real(8), allocatable:: &
       & d_fld(:,:,:), t_fld(:,:,:), ye_fld(:,:,:), en_fld(:,:,:)
  real(8), allocatable:: dx_fld(:,:)
  real(8), dimension(:,:,:,:), allocatable:: x_fld, v_fld, v0_fld

  real(8):: dma(1:npt)

contains


  subroutine set_param(io)

    implicit none

    integer, intent(in):: io


    !..calculation Parametar form './in.dat'
    read(io,*)
    read(io,*) k_zoku, i_test, last_lp, int_t, int_x, nin, nou
    read(io,*)
    read(io,*)
    read(io,*) x1_in, x1_out
    read(io,*)
    read(io,*)
    read(io,*) x2_in, x2_out
    read(io,*)
    read(io,*)
    read(io,*) x3_in, x3_out
    read(io,*)
    read(io,*)
    read(io,*) nout_lpt, n_anim

    close(io)

    !..check
    if( x1_in > x1_out ) stop 'error: bad data #1 set_param @mod_set'
    if( x2_in > x2_out ) stop 'error: bad data #2 set_param @mod_set'
    if( x3_in > x3_out ) stop 'error: bad data #3 set_param @mod_set'

    !..message
    write(*,'(a20,i10)') 'npt = :', npt
    write(*,'("pt move area :")')
    write(*,'(" x1: ", 1pe10.3, 1x, "-->", e10.3)') x1_in, x1_out
    write(*,'(" x2: ", 1pe10.3, 1x, "-->", e10.3)') x2_in, x2_out
    write(*,'(" x3: ", 1pe10.3, 1x, "-->", e10.3)') x3_in, x3_out

    return

  end subroutine set_param



  subroutine set_data

    use mod_cnst, only: pi

    implicit none

    real, allocatable:: x_fld_in(:,:)
    integer:: nx_max, i, j, k, ierr


    !..grid for hydro result
    !..set grid
    read(50) nx1, nx2, nx3

    nx_max = max( max( nx1, nx2 ), nx3 )

    allocate ( x_fld_in(1:ndim,1:nx_max), &
         & x_fld(ndim,nx1,nx2,nx3), dx_fld(ndim,nx_max), &
         & d_fld(nx1,nx2,nx3), t_fld(nx1,nx2,nx3), ye_fld(nx1,nx2,nx3), &
         & en_fld(nx1,nx2,nx3), &
         & v_fld(ndim,nx1,nx2,nx3), v0_fld(ndim,nx1,nx2,nx3), &
         & stat = ierr)

    if(ierr /= 0) stop '### main: error  #3 ###'

    read(50) x_fld_in(1,1:nx1)
    read(50) x_fld_in(2,1:nx2)
    read(50) x_fld_in(3,1:nx3)

    close(50)


    do k = 1, nx3
       do j = 1, nx2
          do i = 1, nx1
             x_fld(1,i,j,k) = dble( x_fld_in(1,i) )
             x_fld(2,i,j,k) = dble( x_fld_in(2,j) )
             x_fld(3,i,j,k) = dble( x_fld_in(3,k) )
          end do
       end do
    end do

    deallocate( x_fld_in )


    dx_fld(1,1) = x_fld(1,1,1,1)
    dx_fld(2,1) = x_fld(2,1,1,1)
    dx_fld(3,1) = x_fld(3,1,1,1)


    if( nx1 >= 2 ) then
       dx_fld(1,2:nx1) = x_fld(1,2:nx1,1,1) - x_fld(1,1:nx1-1,1,1)
    else
       dx_fld(1,1:nx1) = 0.d0
    end if

    if( nx2 >= 2 ) then
       dx_fld(2,2:nx2) = x_fld(2,1,2:nx2,1) - x_fld(2,1,1:nx2-1,1)
    else
       dx_fld(2,1:nx2) = 0.d0
    end if

    if( nx3 >= 2 ) then
       dx_fld(3,2:nx3) = x_fld(3,1,2:nx2,1) - x_fld(3,1,1:nx2-1,1)
    else
       dx_fld(3,1:nx3) = 0.d0
    end if


    !..message
    write(*,'("----- grid information -----")')
    write(*,'("Spherical Coordinate")')
    write(*,'(a20,i5,2(a2,i5))') &
         & 'x1 x x2 x n3 :', nx1, ' x', nx2, ' x', nx3

    write(*,'(a20,1pe10.2,a3,e10.2)') &
         & 'x1-range :'    , x_fld(1,1,1,1), '->',  x_fld(1,nx1,nx2,nx3)
    write(*,'(a20,1pe10.2,a3,e10.2)') &
         & 'x2-range :'    , x_fld(2,1,1,1), '->',  x_fld(2,nx1,nx2,nx3)
    write(*,'(a20,1pe10.2,a3,e10.2)') &
         & 'x3-range :'    , x_fld(3,1,1,1), '->',  x_fld(3,nx1,nx2,nx3)
    write(*,*)


    return

  end subroutine set_data



end module mod_set
