subroutine ofile

  use mod_set, only: k_zoku, set_param

  implicit none

  ! ------------------------------------------------------------------ !
  !     input                                                          !
  ! ------------------------------------------------------------------ !

  open(15, file = './in/lagpat.in', action = 'read')

  call set_param(15)

  !     input                                                          !
  ! ------------------------------------------------------------------ !


  ! ------------------------------------------------------------------ !
  !     output                                                         !
  ! ------------------------------------------------------------------ !

  !..pt all
  open(40, file = './res/init_part.dat', action = 'write')
  open(41, file = './res/fini_part.dat', action = 'write')

  !..position
  open(61, file = './res/hydr/hydr_00001.dat', &
       & form = 'unformatted', action = 'write')

  !..movie
  open(62, file = './res/anim_set.dat', action = 'write')
  open(63, file = './res/anim_ti.dat' , action = 'write')
  open(64, file = './res/anim/anim_00001.dat', action = 'write')

  !..log files
  open(70, file = './res/part_set.log' , action = 'write')
  open(71, file = './res/condition.log', action = 'write')

  !     output                                                         !
  ! ------------------------------------------------------------------ !



  return

end subroutine ofile
