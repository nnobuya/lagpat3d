subroutine sekig_3D(iflag0,t)

  use mod_unit, only: rho_uni, v_uni, tem_uni, r_uni, t_uni
  use mod_cnst, only: ndim
  use mod_set , only: nx1, nx2, nx3, x_fld, dx_fld, &
       & v_fld, v0_fld, d_fld, t_fld, ye_fld, en_fld, &
       & dir_path, eos_name, mass_name
  use mod_data3D, only: nsub_step, njob, njobs, inode, &
       & lv_min, lp_sta, hoge, lv_trc, jproc, kproc, lproc, &
       & js34, je34, j34s, ks34, ke34, k34s, ls34, le34, l34s, &
       & jc, kc, lc, x, y, z, time, qrho, ye, tem, &
       & vlx, vly, vlz, sen, itable,fn

  implicit none

  integer:: omp_get_thread_num, omp_get_max_threads
  integer, save:: j1ma, k1ma, l1ma, lv_max, jr34_min, jr34_max, &
       & kr34_min, kr34_max, lr34_min, lr34_max, job34

  !..local
  integer:: iflag0, iend_thr, imax, imin, irank, ista_thr
  integer:: iwork1, iwork2
  real(8):: t, dlx0

  integer:: i, j, k, l, ii, nn, nphi, nrr, nth, nx_max
  integer:: j1, j2, j3, k1, k2, k3, l1, l2, l3
  integer:: jg, kg, lg, jr, kr, lr, ierr
  integer:: job31, job32, job33
  integer:: job40, job41, job42, job43, job44, job45, job46
  integer:: jr_max, jr31_min, jr31_max, jr32_max, jr32_min

  integer:: lr32_max, lr32_min, lr33_max, lr33_min, lv
  integer:: lvf
  integer:: max_thr, maxs, my_thr
  integer:: ju3, kr31_max, kr31_min, kr33_max, kr33_min

  !---- nishimura subroutine
  real(8), allocatable:: x_fld_in(:,:)

  !..file names
  character*10:: no_index, no_index2




  if (iflag0 /= 0) then

     if( nsub_step == 1 )then

!    allocate(itable(0:jproc-1,0:kproc-1,0:lproc-1))
!    irank = 0 
!    do l = lp_sta, lproc - 1
!     do k = 0, kproc - 1
!      do j = 0, jproc - 1
!       itable(j,k,l) = irank
!       irank = irank + 1
!      end do
!     end do
!    end do

        write(*,'(" Read Sekiguchi 3D data, job number =", i3, &
             & " FMR level for tracer =", i2)') njob,lv_trc

        !write(fn,'(a56,i2.2)')"/misc/work112/sekgchyi/Tracer/SFHo/Dat/grd_SFHo_135-135_",njob
        write(no_index,'(i2.2)') njob
        fn = trim(adjustl(dir_path)) // trim(adjustl(eos_name)) // '/' // &
             & trim(adjustl(mass_name)) // '/Dat/grd_' &
             & // trim(adjustl(eos_name)) // '_' // trim(adjustl(mass_name)) &
             & // '_' // trim(adjustl(no_index))

        open(38,file=fn,status='old')

        read(38,*)dlx0,ju3,lv_max,lvf,nth,nphi,nrr,maxs,                &
             jr31_min,jr31_max,kr31_min,kr31_max,                  &
             jr32_min,jr32_max,lr32_min,lr32_max,                  &
             kr33_min,kr33_max,lr33_min,lr33_max,                  &
             jr34_min,jr34_max,kr34_min,kr34_max,lr34_min,lr34_max,&
             job31,job32,job33,job34,job40,job41,job42,job43,job44,job45,job46
        close(38)

!    allocate(js34(jr34_min:jr34_max),je34(jr34_min:jr34_max),j34s(jr34_min:jr34_max))
!    allocate(ks34(kr34_min:kr34_max),ke34(kr34_min:kr34_max),k34s(kr34_min:kr34_max))
!    allocate(ls34(lr34_min:lr34_max),le34(lr34_min:lr34_max),l34s(lr34_min:lr34_max))

        do jr = jr34_min, jr34_max
           do kr = kr34_min, kr34_max
              do lr = lr34_min, lr34_max
                 i = itable(jr,kr,lr)
        !write(fn,"(a58,i2.2,a1,i4.4)")'/misc/work112/sekgchyi/Tracer/SFHo/Con3D/xyz_SFHo_135-135_',njob,'_',itable(jr,kr,lr)
                 write(no_index ,'(i2.2)') njob
                 write(no_index2,'(i4.4)') itable(jr,kr,lr)
                 fn = trim(adjustl(dir_path)) // trim(adjustl(eos_name)) // '/' // &
                      & trim(adjustl(mass_name)) // '/Con3D/xyz_' // trim(adjustl(eos_name)) &
                      & // '_' // trim(adjustl(mass_name)) // '_' // trim(adjustl(no_index)) // '_' // trim(adjustl(no_index2))

                 open(134+i,file=fn,form='unformatted',status="old")
                 read(134+i)js34(jr),je34(jr),j34s(jr),ks34(kr),ke34(kr),k34s(kr),ls34(lr),le34(lr),l34s(lr)
              end do
           end do
        end do

        j1    = ( jr34_min - jproc/2 )*inode + js34(jr34_min) - 1
        j2    = ( jr34_max - jproc/2 )*inode + je34(jr34_max) - 1
        j3    = j34s(jr34_min)
        j1ma  = ( j2 - j1 )/j3 + 1
        k1    = ( kr34_min - kproc/2 )*inode + ks34(kr34_min) - 1
        k2    = ( kr34_max - kproc/2 )*inode + ke34(kr34_max) - 1
        k3    = k34s(kr34_min)
        k1ma  = ( k2 - k1 )/k3 + 1
        l1    = ( lr34_min - lproc/2 )*inode + ls34(lr34_min) - 1
        l2    = ( lr34_max - lproc/2 )*inode + le34(lr34_max) - 1
        l3    = l34s(lr34_min)
        l1ma  = ( l2 - l1 )/l3 + 1

        if( j1ma /= nx1 .or. k1ma /= nx2 .or. l1ma /= nx3 )then
           write(*,'(" Invalid grid number in sekig 3D =",6i3)')j1ma,nx1,k1ma,nx2,l1ma,nx3
           stop
        endif

!     allocate(jc(j1:j2),kc(k1:k2),lc(l1:l2),time(job34))
!     allocate(x (j1ma,lv_min:lv_max),y(k1ma,lv_min:lv_max),z(l1ma,lv_min:lv_max))
!     allocate(qrho(j1ma,k1ma,l1ma,lv_min:lv_max),VLX(j1ma,k1ma,l1ma,lv_min:lv_max),&
!              VLY (j1ma,k1ma,l1ma,lv_min:lv_max),VLZ(j1ma,k1ma,l1ma,lv_min:lv_max),&
!              tem (j1ma,k1ma,l1ma,lv_min:lv_max),ye (j1ma,k1ma,l1ma,lv_min:lv_max))

        j1 = 1
        do jr = jr34_min, jr34_max
           do j = js34(jr), je34(jr),j34s(jr)
              jg     = ( jr - jproc/2 )*inode + j - 1
              jc(jg) = j1
              do lv = lv_min, lv_max
                 x(j1,lv) = 2.d0**(lv_max-lv)*dlx0*dble(jg)
              end do
              j1     = j1 + 1
           end do
        end do

        k1 = 1
        do kr = kr34_min, kr34_max
           do k = ks34(kr), ke34(kr),k34s(kr)
              kg     = ( kr - kproc/2 )*inode + k - 1 
              kc(kg) = k1
              do lv = lv_min, lv_max
                 y(k1,lv) = 2.d0**(lv_max-lv)*dlx0*dble(kg)
              end do
              k1     = k1 + 1
           end do
        end do

        l1 = 1
        do lr = lr34_min, lr34_max
           do l = ls34(lr), le34(lr),l34s(lr)
              lg     = ( lr - lproc/2 )*inode + l - 1 
              lc(lg) = l1
              do lv = lv_min, lv_max
                 z(l1,lv) = 2.d0**(lv_max-lv)*dlx0*dble(lg)
              enddo
              l1     = l1 + 1
           end do
        end do

        !   endif
        do nn = 1, job34

           !---- read
           imin = itable(jr34_min,kr34_min,lr34_min)
           imax = itable(jr34_max,kr34_max,lr34_max)
           !$omp parallel private(my_thr,max_thr,ista_thr,iend_thr,iwork1,iwork2,&
           !$omp jg,kg,lg,j1,k1,l1,i,jr,kr,lr)
           my_thr   = omp_get_thread_num()
           max_thr  = omp_get_max_threads()
           iwork1   = ( imax - imin + 1 )/max_thr
           iwork2   = mod( imax - imin + 1, max_thr )
           ista_thr = my_thr*iwork1 + imin + min( my_thr, iwork2 )
           iend_thr = ista_thr + iwork1 -1
           if( iwork2 > my_thr ) iend_thr = iend_thr + 1

           do ii = ista_thr,iend_thr
              lr = int( ii / (jproc*kproc) ) + lp_sta
              kr = int(( ii - jproc*kproc*( lr - lp_sta ) )/jproc)
              jr = ii - jproc*kproc*( lr - lp_sta ) - jproc*kr
              i  = itable(jr,kr,lr)

              read(134+i)time(nn)

              do lv = lv_min,lv_max
                 do l = ls34(lr),le34(lr),l34s(lr)
                    do k = ks34(kr),ke34(kr),k34s(kr)
                       do j = js34(jr),je34(jr),j34s(jr)

                          jg = ( jr - jproc/2 )*inode + j - 1
                          kg = ( kr - kproc/2 )*inode + k - 1
                          lg = ( lr - lproc/2 )*inode + l - 1
                          j1 = jc(jg)
                          k1 = kc(kg)
                          l1 = lc(lg)

                          read(134+i) &
                               & qrho(nn,j1,k1,l1,lv),YE (nn,j1,k1,l1,lv),&
                               & TEM(nn,j1,k1,l1,lv), VLX(nn,j1,k1,l1,lv),&
                               & VLY(nn,j1,k1,l1,lv),VLZ(nn,j1,k1,l1,lv)  &
                               , hoge               ,hoge               ,hoge                 &
                               , hoge               ,hoge               ,hoge                 &
                               , sen(nn,j1,k1,l1,lv),hoge               ,hoge                 &
                               , hoge                                                         &
                               , hoge               ,hoge                                     &
                               , hoge               ,hoge

                          !---- cgs unit
                          qrho(nn,j1,k1,l1,lv) = qrho(nn,j1,k1,l1,lv)*rho_uni
                          vlx (nn,j1,k1,l1,lv) = vlx (nn,j1,k1,l1,lv)*v_uni
                          vly (nn,j1,k1,l1,lv) = vly (nn,j1,k1,l1,lv)*v_uni
                          vlz (nn,j1,k1,l1,lv) = vlz (nn,j1,k1,l1,lv)*v_uni
                          tem (nn,j1,k1,l1,lv) = tem (nn,j1,k1,l1,lv)*tem_uni

                       end do
                    end do
                 end do
              end do

           end do
           !$omp end parallel

           write(*,'(" Read 3D sekig data, job =",i3," sub step =",i3," t [ms] =",es10.3)')njob,nn,time(nn)*t_uni*1.d3

        end do

        do jr = jr34_min, jr34_max
           do kr = kr34_min, kr34_max
              do lr = lr34_min, lr34_max
                 i = itable(jr,kr,lr)
                 close(134+i)
              end do
           end do
        end do

!       deallocate(itable,jc,kc,lc,time,x,y,z)!,qrho,vlx,vly,vlz,tem,ye)
!       deallocate(js34,je34,j34s,ks34,ke34,k34s,ls34,le34,l34s)

     endif

     nsub_step = job34 - nsub_step + 1
     t         = time(nsub_step)*t_uni

     if( nsub_step <= 1 )then

        nsub_step = 1
        njob      = njob - 1
        if( njob < njobs )then
           write(*,'(" No more fluid data")')
           stop
        endif
 
     endif

  endif

   !---- from nishimura set_data (excute only for initial)
  if( iflag0 == 0 )then

     allocate(itable(0:jproc-1,0:kproc-1,0:lproc-1))
     irank = 0 
     do l = lp_sta, lproc - 1
        do k = 0, kproc - 1
           do j = 0, jproc - 1
              itable(j,k,l) = irank
              irank = irank + 1
           end do
        end do
     end do

     write(*,'(" Read Sekiguchi 3D data, FMR level for tracer =",i2)')lv_trc

     !write(fn,'(a56,i2.2)')"/misc/work112/sekgchyi/Tracer/SFHo/Dat/grd_SFHo_135-135_",njob

     write(no_index,'(i2.2)') njob
     fn = trim(adjustl(dir_path)) // trim(adjustl(eos_name)) // '/' // &
          & trim(adjustl(mass_name)) // '/Dat/grd_' // trim(adjustl(eos_name)) &
          & // '_' // trim(adjustl(mass_name)) // '_' // trim(adjustl(no_index))

     open(38,file=fn,status='old')
     read(38,*)dlx0,ju3,lv_max,lvf,nth,nphi,nrr,maxs,                &
          jr31_min,jr31_max,kr31_min,kr31_max,                  &
          jr32_min,jr32_max,lr32_min,lr32_max,                  &
          kr33_min,kr33_max,lr33_min,lr33_max,                  &
          jr34_min,jr34_max,kr34_min,kr34_max,lr34_min,lr34_max,&
          job31,job32,job33,job34,job40,job41,job42,job43,job44,job45,job46
     close(38)

     job34 = 10

     allocate(js34(jr34_min:jr34_max),je34(jr34_min:jr34_max),j34s(jr34_min:jr34_max))
     allocate(ks34(kr34_min:kr34_max),ke34(kr34_min:kr34_max),k34s(kr34_min:kr34_max))
     allocate(ls34(lr34_min:lr34_max),le34(lr34_min:lr34_max),l34s(lr34_min:lr34_max))

     do jr = jr34_min, jr34_max
        do kr = kr34_min, kr34_max
           do lr = lr34_min, lr34_max
              i = itable(jr,kr,lr)
              !write(fn,"(a58,i2.2,a1,i4.4)")'/misc/work112/sekgchyi/Tracer/SFHo/Con3D/xyz_SFHo_135-135_',njob,'_',itable(jr,kr,lr)

              write(no_index ,'(i2.2)') njob
              write(no_index2,'(i4.4)') itable(jr,kr,lr)
              fn = trim(adjustl(dir_path)) // trim(adjustl(eos_name)) // '/' // &
                   & trim(adjustl(mass_name)) // '/Con3D/xyz_' // trim(adjustl(eos_name)) &
                   & // '_' // trim(adjustl(mass_name)) // '_' // trim(adjustl(no_index)) // '_' // trim(adjustl(no_index2))

              open(134+i,file=fn,form='unformatted',status="old")
              read(134+i)js34(jr),je34(jr),j34s(jr),ks34(kr),ke34(kr),k34s(kr),ls34(lr),le34(lr),l34s(lr)
           end do
        end do
     end do

     j1    = ( jr34_min - jproc/2 )*inode + js34(jr34_min) - 1
     j2    = ( jr34_max - jproc/2 )*inode + je34(jr34_max) - 1
     j3    = j34s(jr34_min)
     j1ma  = ( j2 - j1 )/j3 + 1
     k1    = ( kr34_min - kproc/2 )*inode + ks34(kr34_min) - 1
     k2    = ( kr34_max - kproc/2 )*inode + ke34(kr34_max) - 1
     k3    = k34s(kr34_min)
     k1ma  = ( k2 - k1 )/k3 + 1
     l1    = ( lr34_min - lproc/2 )*inode + ls34(lr34_min) - 1
     l2    = ( lr34_max - lproc/2 )*inode + le34(lr34_max) - 1
     l3    = l34s(lr34_min)
     l1ma  = ( l2 - l1 )/l3 + 1

     allocate(jc(j1:j2),kc(k1:k2),lc(l1:l2))
     allocate(x (j1ma,lv_min:lv_max),y(k1ma,lv_min:lv_max),z(l1ma,lv_min:lv_max))
     allocate(time(job34))
     allocate(qrho(job34,j1ma,k1ma,l1ma,lv_min:lv_max),VLX(job34,j1ma,k1ma,l1ma,lv_min:lv_max),&
          VLY (job34,j1ma,k1ma,l1ma,lv_min:lv_max),VLZ(job34,j1ma,k1ma,l1ma,lv_min:lv_max),&
          tem (job34,j1ma,k1ma,l1ma,lv_min:lv_max),ye (job34,j1ma,k1ma,l1ma,lv_min:lv_max),&
          sen (job34,j1ma,k1ma,l1ma,lv_min:lv_max))

     j1 = 1
     do jr = jr34_min, jr34_max
        do j = js34(jr), je34(jr),j34s(jr)
           jg     = ( jr - jproc/2 )*inode + j - 1
           jc(jg) = j1
           do lv = lv_min, lv_max
              x(j1,lv) = 2.d0**(lv_max-lv)*dlx0*dble(jg)
           end do
           j1     = j1 + 1
        end do
     end do

     k1 = 1
     do kr = kr34_min, kr34_max
        do k = ks34(kr), ke34(kr),k34s(kr)
           kg     = ( kr - kproc/2 )*inode + k - 1 
           kc(kg) = k1
           do lv = lv_min, lv_max
              y(k1,lv) = 2.d0**(lv_max-lv)*dlx0*dble(kg)
           end do
           k1     = k1 + 1
        end do
     end do

     l1 = 1
     do lr = lr34_min, lr34_max
        do l = ls34(lr), le34(lr),l34s(lr)
           lg     = ( lr - lproc/2 )*inode + l - 1 
           lc(lg) = l1
           do lv = lv_min, lv_max
              z(l1,lv) = 2.d0**(lv_max-lv)*dlx0*dble(lg)
           enddo
           l1     = l1 + 1
        end do
     end do

     do jr = jr34_min, jr34_max
        do kr = kr34_min, kr34_max
           do lr = lr34_min, lr34_max
              i = itable(jr,kr,lr)
              close(134+i)
           end do
        end do
     end do

     !..grid for hydro result
     !..set grid
     nx1 = j1ma
     nx2 = k1ma
     nx3 = l1ma
!    read(50) nx1, nx2, nx3

     nx_max = max( max( nx1, nx2 ), nx3 )
     allocate ( x_fld_in(1:ndim,1:nx_max), &
          & x_fld(ndim,nx1,nx2,nx3), dx_fld(ndim,nx_max), &
          & d_fld(nx1,nx2,nx3), t_fld(nx1,nx2,nx3), ye_fld(nx1,nx2,nx3), &
          & en_fld(nx1,nx2,nx3), &
          & v_fld(ndim,nx1,nx2,nx3), v0_fld(ndim,nx1,nx2,nx3), &
          & stat = ierr)

     if( ierr /= 0 ) stop '### main: error  #3 ###'

     do j = 1, nx1
        x_fld_in(1,j) = x(j,lv_trc)*r_uni
     end do

     do k = 1, nx2
        x_fld_in(2,k) = y(k,lv_trc)*r_uni
     end do

     do l = 1, nx3
        x_fld_in(3,l) = z(l,lv_trc)*r_uni
     end do

!    read(50) x_fld_in(1,1:nx1)
!    read(50) x_fld_in(2,1:nx2)
!    read(50) x_fld_in(3,1:nx3)

!    close(50)

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

!    if ( i_test == 1 ) then
!       write(80,*) nx1, nx2, nx3
!       do k = 1, nx3
!          do j = 1, nx2
!             do i = 1, nx1
!                write(80,*) x_fld(1,i,j,k), x_fld(2,i,j,k), x_fld(3,i,j,k)
!             end do
!          end do
!       end do
!    end if


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


!    deallocate(itable,jc,kc,lc,x,y,z)
!    deallocate(js34,je34,j34s,ks34,ke34,k34s,ls34,le34,l34s)

  end if


  return

end subroutine sekig_3D
