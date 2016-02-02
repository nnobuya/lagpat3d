module mod_data3d

  implicit none

  private
  public:: nsub_step, njob, &
       & lv_min, lp_sta, hoge, lv_trc, jproc, kproc, lproc, &
       & js34, je34, j34s, ks34, ke34, k34s, ls34, le34, l34s, &
       & jc, kc, lc, x, y, z, time, qrho, ye, tem, &
       & vlx, vly, vlz, sen, itable,fn

  integer njob, nsub_step
  integer,allocatable :: itable(:,:,:)
  character(100):: fn
  integer, allocatable:: js34(:),je34(:),j34s(:)
  integer, allocatable:: ks34(:),ke34(:),k34s(:)
  integer, allocatable:: ls34(:),le34(:),l34s(:)
  integer, allocatable:: jc(:),kc(:),lc(:)
  real(4), allocatable:: x(:,:),y(:,:),z(:,:),time(:)
  real(4), allocatable:: qrho(:,:,:,:,:),VLX(:,:,:,:,:),VLY(:,:,:,:,:), &
  & VLZ(:,:,:,:,:), ye(:,:,:,:,:), tem(:,:,:,:,:), sen(:,:,:,:,:)
  real(4) hoge

!  !--- xyz SFHo 1.35-1.35
!  integer, parameter:: njobs  = 10, njobe = 99, inode  = 78, &
!       & jproc  = 4, kproc  = jproc, lproc  = jproc, lp_sta = lproc/2

  !--- based on "xyz SFHo 1.35-1.35"
  integer, parameter:: &
       & jproc  = 4, kproc  = jproc, lproc  = jproc, lp_sta = lproc/2
  integer, parameter:: lv_min = 1 , lv_trc = 2


end module mod_data3d
