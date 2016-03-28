program traj_rev

  implicit none

  !..param.
  integer  , parameter:: ndim = 3
  character, parameter:: bin*11 = 'unformatted'
  real(4)  , parameter:: r_mev  = 1.16e10
  real(4)  , parameter:: te_nse = 0.5 *r_mev

  integer, parameter:: npt0  = 100
  logical, parameter:: debug = .true.

  integer:: npt, idt_ini, idt_fin

  integer, allocatable:: istg(:), ipt(:,:,:), ist_pt(:,:)
  real(4), allocatable:: ti(:), dt(:), x_pt(:,:,:), v_pt(:,:,:), &
      & d_pt(:,:), t_pt(:,:), ye_pt(:,:), en_pt(:,:)
  real(4):: rad

  integer:: idt, ier, i, j, k
  character:: op_file*100


  !..open files
  open(50,file = './res/part_set.log', action = 'read')


  read(50,*)
  read(50,*) npt, i, j, k
  close(50)

  idt_ini = 1
  idt_fin = 162

  allocate( istg(idt_ini:idt_fin), ti(idt_ini:idt_fin), dt(idt_ini:idt_fin), &
       & ist_pt(1:npt,idt_ini:idt_fin), ipt(1:ndim,1:npt,idt_ini:idt_fin), &
       & x_pt(1:ndim,1:npt,idt_ini:idt_fin), &
       & v_pt(1:ndim,1:npt,idt_ini:idt_fin), d_pt(1:npt,idt_ini:idt_fin), &
       & t_pt(1:npt,idt_ini:idt_fin), ye_pt(1:npt,idt_ini:idt_fin), &
       & en_pt(1:npt,idt_ini:idt_fin), stat = ier)

  write(*,'(5x, "- read data")')

  lp_read_hydr: do idt = idt_ini, idt_fin

     if (mod(idt,50) == 0) write(*,*) idt

     write(op_file,'("./res/hydr/hydr_", i5.5, ".dat")') idt

     open(50,file = op_file, form = bin, action = 'read')

     read(50) istg(idt), ti(idt), dt(idt)
     read(50) ist_pt(1:npt,idt)
     read(50) ipt(1:ndim,1:npt,idt), x_pt(1:ndim,1:npt,idt), &
          & v_pt(1:ndim,1:npt,idt), d_pt(1:npt,idt), t_pt(1:npt,idt), &
          & ye_pt(1:npt,idt), en_pt(1:npt,idt)

     close(50)

  end do lp_read_hydr



  write(*,'(5x, "- write data")')

  if (debug) npt = npt0

  lp_write_data: do i = 1, npt

     if (mod(i,5000) == 0) write(*,*) i, npt

     write(op_file,'("./res/traj/traj_", i7.7, ".dat")') i
     open(60,file = op_file, action = 'write')

     write(60,'("#   Time", 11x, "Density", 8x, "T", 14x, &
          & "Entropy", 8x, "Ye", 13x, "Radius     ")')
     lp_pt_evol:do idt = idt_fin, idt_ini, -1
        rad = sqrt(sum(x_pt(1:ndim,i,idt) *x_pt(1:ndim,i,idt)))
        write(60,'(1p, *(e15.7))') &
           & ti(idt) - ti(idt_fin), &
           & d_pt(i,idt), t_pt(i,idt), en_pt(i,idt), ye_pt(i,idt), rad
     end do lp_pt_evol
     close(60)


     !..Ye-S map at 0.5 MeV

     if (maxval(t_pt(i,idt_ini:idt_fin)) < te_nse) then
        print *, t_pt(i,idt_fin), maxval(t_pt(i,idt_ini:idt_fin))
     else
        lp_search: do idt = idt_fin, idt_ini, -1
           if (t_pt(i,idt) >= te_nse) exit lp_search
        end do lp_search

        print *, idt

     end if





  end do lp_write_data





  stop 'traj_rev: normal termination'

end program traj_rev
