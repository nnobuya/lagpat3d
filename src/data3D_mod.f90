module data3D
 implicit real(8) (a-h,o-z)
 integer njob,nsub_step
 integer,allocatable :: itable(:,:,:)
 character(100) fn
 integer,allocatable :: js34(:),je34(:),j34s(:)
 integer,allocatable :: ks34(:),ke34(:),k34s(:)
 integer,allocatable :: ls34(:),le34(:),l34s(:)
 integer,allocatable :: jc(:),kc(:),lc(:)
 real(4),allocatable :: x(:,:),y(:,:),z(:,:),time(:)
 real(4),allocatable :: qrho(:,:,:,:,:),VLX(:,:,:,:,:),VLY(:,:,:,:,:),VLZ(:,:,:,:,:),&
                        ye  (:,:,:,:,:),tem(:,:,:,:,:),sen(:,:,:,:,:)
 real(4) hoge
 !--- xyz SFHo 1.35-1.35
 parameter( njobs  = 10, njobe = 99, inode  = 78, jproc  = 4, kproc  = jproc, lproc  = jproc, lp_sta = lproc/2 )
 parameter( lv_min = 1 , lv_trc = 2 )
end module data3D
