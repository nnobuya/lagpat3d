program traj_rev

  implicit none

  integer  , parameter:: ndim = 3
  character, parameter:: bin*11 = 'unformatted'

  integer:: npt, idt_ini, idt_fin
  integer, allocatable::

  integer:: idt, i, j, k
  character:: op_file*100


  open(50,file = './res/part_set.log', action = 'read')
  read(50,*)
  read(50,*) npt, i, j, k
  close(50)

  idt_ini = 1
  idt_fin = 162


  lp_read_hydr: do idt = idt_ini, idt_fin

     write(op_file,'("./res/hydr/hydr_", i5.5, ".dat")') idt

     open(50,file = op_file, form = bin, action = 'read')

     read(50) istg(idt), ti(idt), dt(idt)
     read(50) ist_pt(1:npt,idt)
     read(50) ipt(1:ndim,1:npt), x_pt(1:ndim,1:npt), v_pt(1:ndim,1:npt), &
          & d_pt(1:npt), t_pt(1:npt), ye_pt(1:npt), en_pt(1:npt)


     close(50)

  end do lp_read_hydr



  stop 'traj_rev: normal termination'

end program traj_rev
