! -------------------------------------------------------------------- !
!                                                                      !
! Lag_rangian Pa_rticle T_racer for 3-D                                !
!                                                                      !
! Purpose: to trace motions of particles from hydro results            !
! time stamp: 2005/09/24                                               !
!             2006/01/04                                               !
!             2009/04/16                                               !
!             2009/08/25                                               !
!             2010/01/24                                               !
!             2012/07/24                                               !
!             2013/01/28                                               !
!             2014/10/21                                               !
!                                                                      !
! input : hydro Result: LAGPAT/lag_*.dat => Sekiguchi 3D data          !
!         ./in.dat                                                     !
!                                                                      !
! output: ./part/*.dat                                                 !
!                                                                      !
! Coded by Nishimura Nobuya                                            !
!                                                                      !
! -------------------------------------------------------------------- !

program lagpat

  use mod_cnst, only: npt, ndim, set_cnst
  use mod_set , only: int_t, i_test, last_lp, &
       & d_fld, t_fld, ye_fld, en_fld, v_fld, v0_fld, &
       & set_data
  use mod_fld   , only: dt_max, fld
  use mod_data3D, only: njob, njobe, nsub_step

  implicit none

  !..main
  integer, allocatable:: ist_pt(:), ipt(:,:), id(:,:)
  real(8), allocatable:: d_pt(:), t_pt(:), ye_pt(:), en_pt(:), &
       & x_pt(:,:), v_pt(:,:), v_pt_p(:,:,:), dma(:)

  !..local
  integer:: istg, n_anim_out = 0
  integer:: ier
  real(8):: ti, dt0, dt_in, dt = 0.d0


  ! ------------------------------------------------------------------ !
  !                                                                    !
  !     pre-process                                                    !
  !                                                                    !
  ! ------------------------------------------------------------------ !

  !..logo
  call logo

  !..open files
  call ofile


  !..allocation

  call set_cnst

  allocate(ist_pt(1:npt), ipt(1:ndim,1:npt), id(1:ndim,1:npt), &
       & d_pt(1:npt), t_pt(1:npt), ye_pt(1:npt), en_pt(1:npt), &
       & x_pt(1:ndim,1:npt), v_pt(1:ndim,1:npt), &
       & v_pt_p(1:ndim,0:4,1:npt), dma(1:npt), stat=ier)
  if (ier /= 0) stop 'lagpat(): allocation error'


  !..set field data
  njob      = njobe
  nsub_step = 1
  call sekig_3D(0,ti)


  close(41)
  close(42)
  close(60)
  close(91)


  call fld(ier, d_fld(:,:,:), t_fld(:,:,:), ye_fld(:,:,:), en_fld(:,:,:), &
       & v_fld(:,:,:,:), v0_fld(:,:,:,:))
  ! out: all

  call set_pt(istg, ti, ist_pt(:), id(:,:), dma(:), &
  & x_pt(:,:), v_pt(:,:), d_fld(:,:,:))
  !  out: all

  call hokan_main(1, dt_max, ist_pt(:), ipt(:,:), x_pt(:,:), &
       & d_fld(:,:,:), t_fld(:,:,:), ye_fld(:,:,:), en_fld(:,:,:), &
       & v0_fld(:,:,:,:), v_fld(:,:,:,:), &
       & d_pt(:), t_pt(:), ye_pt(:), en_pt(:), v_pt(:,:), v_pt_p(:,:,:))

  !  in: mode, ist_pt, rad_pt, the_pt, d_fld, t_fld, v_fld
  !  in: d_pt, t_pt, v_pt
  !  inout: mode, ipt, jpt

  dt0 = dt_max
  dt  = dt_max

  if( i_test == 1 ) stop '### finish test  ###'


  !     pre-process                                                    !
  ! ------------------------------------------------------------------ !


  ! ------------------------------------------------------------------ !
  !                                                                    !
  !     particle tracing: LagPaT main                                  !
  !                                                                    !
  ! ------------------------------------------------------------------ !

  main_lp: do

     istg = istg + 1

     !..sepcial for steady flow
     !dt0 = 1.d-4
     !dt  = 1.d-4
     !..end


     ! --------------------------------------------------------------- !
     !     output                                                      !
     ! --------------------------------------------------------------- !

     call output( istg, ti, dt, ipt(:,:), ist_pt(:), &
          & x_pt(:,:), d_pt(:), t_pt(:), ye_pt(:), en_pt(:), v_pt(:,:), &
          & n_anim_out )
     !    in: others
     ! inout: n_anim_out

     !     output                                                      !
     ! --------------------------------------------------------------- !


     ! --------------------------------------------------------------- !
     !     update                                                      !
     ! --------------------------------------------------------------- !


!NN     call dt_chk( ti, dt, ist_pt(:), ipt(:,:), v_pt(:,:), x_pt(:,:) )
     !  in: all

     call move( dt, dt0, ist_pt(:), v_pt(:,:), v_pt_p(:,:,:), x_pt(:,:) )
     !  in : dt, ist_pt, ipt, jpt, v_pt
     !  out: rad_pt, the_pt

     call status( x_pt(:,:), ist_pt(:) )
     !  in: rad_pt
     ! i&o: ist_pt


     !     update                                                      !
     ! --------------------------------------------------------------- !


     ! --------------------------------------------------------------- !
     !     read                                                        !
     ! --------------------------------------------------------------- !

     call fld(ier, d_fld(:,:,:), t_fld(:,:,:), ye_fld(:,:,:), &
          & en_fld(:,:,:), v_fld(:,:,:,:), v0_fld(:,:,:,:))
     !  out: all

     if( ier /= 0 .or. (last_lp > 0 .and. istg > last_lp) ) exit main_lp

     !     read                                                        !
     ! --------------------------------------------------------------- !

     dt0 = dt
     dt  = dt_max

     ti = ti + dt

     ! --------------------------------------------------------------- !
     !     hokan                                                       !
     ! --------------------------------------------------------------- !

     if ( int_t == 2 ) then
        dt_in = dt0
     else
        dt_in = dt
     end if

     call hokan_main( 2, dt_in, ist_pt(:), ipt(:,:), x_pt(:,:), &
          & d_fld(:,:,:), t_fld(:,:,:), ye_fld(:,:,:), en_fld(:,:,:), &
          & v0_fld(:,:,:,:), v_fld(:,:,:,:), &
          & d_pt(:), t_pt(:), ye_pt(:), en_pt(:), v_pt(:,:), v_pt_p(:,:,:) )

     !    in: others
     !   out: d_pt, t_pt, v_pt
     ! inout: ipt, jpt


     !     hokan                                                       !
     ! --------------------------------------------------------------- !
     
  end do main_lp

  close(51)
  close(52)
  close(53)
  close(54)
  close(55)
  close(56)
  close(57)

  close(61)
  close(62)
  close(63)

  !     particle tracing: LagPaT main                                  !
  ! ------------------------------------------------------------------ !



  ! ------------------------------------------------------------------ !
  !     closing                                                        !
  ! ------------------------------------------------------------------ !

  write(65,'(a20,i10)') 'calculation step:', istg
  write(65,'(a20,i10)') 'output:', n_anim_out
  close(65)

  write(*,'(a20,i10)') 'calculation step:', istg
  write(*,'(a20,i10)') 'output:', n_anim_out

  call fini_out(istg, ti, ist_pt(:), id(:,:), dma(:), x_pt(:,:), v_pt(:,:))

  close(90)


  !     closing                                                        !
  ! ------------------------------------------------------------------ !


  stop 'lagpat: Normal Termination'


end program lagpat
