subroutine output(istg, ti, dt, ipt, &
     & ist_pt, ist_at, x_pt, d_pt, t_pt, ye_pt, en_pt, v_pt, n_anim_out )

  use mod_cnst, only: npt, ndim
  use mod_set , only: nout_lpt, n_anim

  implicit none

  !..io
  integer, intent(in):: istg, ist_pt(1:npt), ist_at(1:npt), ipt(1:ndim,1:npt)
  real(8), intent(in):: ti, dt, x_pt(1:ndim,npt), v_pt(ndim,npt), &
       & d_pt(npt), t_pt(npt), ye_pt(npt), en_pt(npt)
  integer, intent(inout):: n_anim_out

  !..local
  integer:: npt_in, npt_out, npt_num
  integer, save:: i_hydr = 1
  integer, save:: i_anim = 1
  integer:: i, jpt
  character:: op_file*100


  !..physical quantities
  out_cond: if( istg == 1 .or. mod(istg,nout_lpt) == 0 ) then
     write(op_file,'("../res/hydr/hydr.", i7.7, ".dat")') i_hydr
     open(61, file = op_file, form = 'unformatted', action = 'write')
     !..write
     write(61) istg, real(ti), real(dt)
     write(61) ist_pt(1:npt)
     write(61) ipt(1:ndim,1:npt), &
          & real(x_pt(1:ndim,1:npt)), &
          & real(v_pt(1:ndim,1:npt)), &
          & real( d_pt(1:npt)), real( t_pt(1:npt)), &
          & real(ye_pt(1:npt)), real(en_pt(1:npt))
     close(61)
     i_hydr = i_hydr + 1
  end if out_cond


  if ( mod(istg,n_anim) == 0 ) then
     write(63,'(f10.5)') ti *1.d3

     write(op_file,'("../res/anim/anim.", i7.7, ".dat")') i_anim
     !..write
     open(64, file = op_file, action = 'write')
     do jpt = 1, npt
        if ( ist_pt(jpt) /= 0 ) cycle
        write(64,'(3f9.2,f8.4,i10)') x_pt(1:ndim,jpt) *1.d-5, ye_pt(jpt), jpt
     end do
     close(64)
     i_anim = i_anim + 1

     n_anim_out = n_anim_out + 1
  end if

  call pt_status( ist_pt(:), npt_in, npt_out, npt_num )

  if( istg == 1 ) then
     write(*,*) &
          & '======================= tracing particles  ======================='
     write(*,'("#", 5x, "step", 8x, "ti", 10x, "dt", &
          & 3x, "pt(move)", 4x, "pt(in)", 3x, "pt(out)")')
  end if

  if( mod(istg,10) == 0 .or. istg <= 5 ) &
       & write(*,'(i10,1p2e12.4,3i10)') &
       & istg, ti, dt, npt_num, npt_in, npt_out
!  if( mod(istg,10000) == 0 .or. istg <= 5 ) &
!       & write(*,'(i10,1p2e12.4,3i10)') &
!       & istg, ti, dt, npt_num, npt_in, npt_out

  return

end subroutine output
