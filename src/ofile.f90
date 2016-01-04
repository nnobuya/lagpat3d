subroutine ofile

  use mod_set, only: k_zoku, set_param

  implicit none

  character:: path*50

  character*60:: grid, v1, v2, v3, de, te, ye, en

  ! ------------------------------------------------------------------ !
  !     input                                                          !
  ! ------------------------------------------------------------------ !

  open(15, file = '../in/lagpat.in', action = 'read')

  read(15,*)
  read(15,*) path
  read(15,*)

  call set_param(15)

  !     input                                                          !
  ! ------------------------------------------------------------------ !


  ! ------------------------------------------------------------------ !
  !     output                                                         !
  ! ------------------------------------------------------------------ !

  !..settings
  open(41, file = '../res/part_mass.log', action='write')
  open(42, file = '../res/init_part.dat', action='write')

  !..particle motion
  open(60, file = '../res/set.dat'  , action = 'write')

  !! position
  open(61, file = '../res/hydr/hydr_lpt.0000001.dat', &
       & form = 'unformatted', action = 'write')

  !..movie
  open(62, file = '../res/anim_set.dat', action = 'write')
  open(63, file = '../res/anim_ti.dat' , action = 'write')
  open(64, file = '../res/anim/anim.0000001.dat', action = 'write')


  !..log files
  open(70, file = '../res/condition.log', action = 'write')


  !! final status
  open(90, file = '../res/fini.dat', action = 'write')

  open(92, file = '../res/set.ns.dat' , action = 'write')
  open(93, file = '../res/fini.ns.dat', action = 'write')


  !     output                                                         !
  ! ------------------------------------------------------------------ !



  return

end subroutine ofile
