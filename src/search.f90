subroutine search(mode, xpt, vpt, fac, ipt)

  use mod_cnst, only: ndim
  use mod_set , only: nx1, nx2, nx3, x_fld

  implicit none

  !..main
  integer         , intent(in)   :: mode
  real(8), intent(in)   :: xpt(1:ndim), vpt(1:ndim)
  integer         , intent(inout):: ipt(1:ndim)
  real(8), intent(out)  :: fac(1:ndim)

  !..temp
  integer:: nd
  integer:: ipt_in(1:ndim)
  integer:: i, j, k, ier


  ipt_in(1:ndim) = ipt(1:ndim)


  if( mode == 1 ) then

     if( nx1 >= 2 ) then
        ipt(1) = 1
        j      = 1
        k      = 1
        do i = 1, nx1 - 1
           if( x_fld(1,i,j,k) >= xpt(1) ) exit
           ipt(1) = i
        end do
        fac(1) = ( xpt(1) - x_fld(1,ipt(1),j,k) ) &
             & /( x_fld(1,ipt(1)+1,j,k) - x_fld(1,ipt(1),j,k) )
     else
        ipt(1) = 1
        fac(1) = 0.d0
     end if


     if( nx2 >= 2 ) then
        i      = 1
        ipt(2) = 1
        k      = 1
        do j = 1, nx2 - 1
           if( x_fld(2,i,j,k) >= xpt(2) ) exit
           ipt(2) = j
        end do
        fac(2) = ( xpt(2) - x_fld(2,i,ipt(2),k) ) &
             & /( x_fld(2,i,ipt(2)+1,k) - x_fld(2,i,ipt(2),k) )
     else
        ipt(2) = 1
        fac(2) = 0.d0
     end if

     if( nx3 >= 2 ) then
        i      = 1
        j      = 1
        ipt(3) = 1
        do k = 1, nx3 - 1
           if ( x_fld(3,i,j,k) >= xpt(3) ) exit
           ipt(3) = k
        end do
        fac(3) = ( xpt(3) - x_fld(3,i,j,ipt(3)) ) &
             & /( x_fld(3,i,j,ipt(3)+1) - x_fld(2,i,j,ipt(3)) )
     else
        ipt(3) = 1
        fac(3) = 0.d0
     end if


  else if( mode == 2 ) then

     if( nx1 >= 2 ) then
        j = 1
        k = 1
        if ( xpt(1) >= x_fld(1,nx1 - 1,j,k) ) then
           ipt(1) = nx1 - 1
        else if ( xpt(1) < x_fld(1,1,j,k) ) then
           ipt(1) = 1
        else
           if ( vpt(1) >= 0 ) nd =  1
           if ( vpt(1) <  0 ) nd = -1
           ier = 0
           lp_search_i1: do
              if ( ( ipt(1) == 0 ) .or. ( ipt(1) == nx1 ) ) then
                 ier = - 1
                 exit lp_search_i1
              else if( x_fld(1,ipt(1),j,k) <= xpt(1) &
                   &  .and. xpt(1) < x_fld(1,ipt(1) + 1,j,k) ) then
                 exit lp_search_i1
              end if

              ipt(1) = ipt(1) + nd

           end do lp_search_i1


           if ( ier /= 0 ) then
!NN              write(*,*) 'warning: search radius-dir'
              ipt(1) = ipt_in(1)
              ier = 0
              lp_search_i12: do
                 if( x_fld(1,ipt(1),j,k) <= xpt(1) &
                      & .and. xpt(1) < x_fld(1,ipt(1) + 1,j,k) ) exit lp_search_i12
                 ipt(1) = ipt(1) - nd
                 if ( ( ipt(1) == 0 ) .or. ( ipt(1) == nx1 ) ) ier = - 1
              end do lp_search_i12
           end if

           if ( ier /= 0 ) stop 'search(): error #1'

        end if

        ipt(1) = min( max(ipt(1),1), nx1 - 1 )
        fac(1) = ( xpt(1) - x_fld(1,ipt(1),j,k) ) &
             & /( x_fld(1,ipt(1) + 1,j,k) - x_fld(1,ipt(1),j,k) )
     else
        ipt(1) = 1
        fac(1) = 0.d0
     end if


     if( nx2 >= 2 ) then
        i = 1
        k = 1
        if ( xpt(2) >= x_fld(2,i,nx2-1,k)  ) then
           ipt(2) = nx2 - 1
        else if ( xpt(2) < x_fld(2,i,1,k) ) then
           ipt(2) = 1
        else
           if ( vpt(2) >= 0 ) nd =  1
           if ( vpt(2) <  0 ) nd = -1
           ier = 0
           lp_search_i2: do
              if ( ( ipt(2) == 0 ) .or. ( ipt(2) == nx2 ) ) then
                 ier = - 1
                 exit lp_search_i2
              else if ( ( x_fld(2,i,ipt(2),k) <= xpt(2) ) &
                   & .and. ( xpt(2) < x_fld(2,i,ipt(2)+1,k) ) ) then
                 exit lp_search_i2
              end if
              ipt(2) = ipt(2) + nd
           end do lp_search_i2

           if ( ier /= 0 ) then

              !NN write(*,*) 'warning: search theta-dir'

              ipt(2) = ipt_in(2)
              ier = 0
              lp_search_i22: do

                 if ( ( ipt(2) == 0 ) .or. ( ipt(2) == nx2 ) ) then
                    ier = - 1
                    exit lp_search_i22
                 else if ( ( x_fld(2,i,ipt(2),k) <= xpt(2) ) &
                      &  .and. ( xpt(2) < x_fld(2,i,ipt(2)+1,k) ) ) then
                    exit lp_search_i22
                 end if
                 ipt(2) = ipt(2) - nd
              end do lp_search_i22
           end if

           if ( ier /= 0 ) stop 'search(): error #2'

        end if

        ipt(2) = min( max(ipt(2),1), nx2 - 1 )

        fac(2) = ( xpt(2) - x_fld(2,i,ipt(2),k) ) &
             & /( x_fld(2,i,ipt(2)+1,k) - x_fld(2,i,ipt(2),k) )
     else
        ipt(2) = 1
        fac(2) = 0.d0
     end if


     if ( nx3 >= 2 ) then
        i = 1
        j = 1
        if ( xpt(3) >= x_fld(3,i,j,nx3 - 1)  ) then
           ipt(3) = nx3 - 1
        else if ( xpt(3) <= x_fld(3,i,j,1) ) then
           ipt(3) = 1
        else
           if ( vpt(3) >= 0 ) nd =  1
           if ( vpt(3) <  0 ) nd = -1
           ier = 0
           lp_search_i3: do
              if ( ( ipt(3) == 0 ) .or. ( ipt(3) == nx3 ) ) then
                 ier = - 1
                 exit lp_search_i3
              else if ( ( x_fld(3,i,j,ipt(3)) <= xpt(3) ) &
                   & .and. ( xpt(3) < x_fld(2,i,j,ipt(3)+1) ) ) then
                 exit lp_search_i3
              end if
              ipt(3) = ipt(3) + nd
           end do lp_search_i3

           if ( ier /= 0 ) then

              !NN write(*,*) 'warning: search x3-dir'

              ipt(3) = ipt_in(3)
              ier = 0
              lp_search_i32: do
                 if ( ( ipt(3) == 0 ) .or. ( ipt(3) == nx3 ) ) then
                    ier = - 1
                    exit lp_search_i32
                 else if ( ( x_fld(3,i,j,ipt(3)) <= xpt(3) ) &
                      &  .and. ( xpt(3) < x_fld(3,i,j,ipt(3)+1) ) ) then
                    exit lp_search_i32
                 end if
                 ipt(3) = ipt(3) - nd
              end do lp_search_i32
           end if

           if ( ier /= 0 ) stop 'search(): error #3'
!NN 2014/10/23           if ( ier /= 0 ) stop 'search(): error #3'

        end if

        ipt(3) = min( max(ipt(3),1), nx3 - 1 )
        fac(3) = ( xpt(3) - x_fld(3,i,j,ipt(3)) ) &
             & /( x_fld(3,i,j,ipt(3)+1) - x_fld(2,i,j,ipt(3)) )

     else
        ipt(3) = 1
        fac(3) = 0.d0
     end if

  else
     stop 'search(): error'
  end if


  return

end subroutine search
