! One name for the two tables below, so `cint_stg_roots` needs only this
! module and does not care that the data is split across two files.
!
! The split is GitHub's 100 MB per-file limit and nothing deeper: the two
! arrays together are 108 MB of source, and separately they are 54 MB each.
module cint_tab_stg_roots
   use cint_tab_stg_x, only: STG_X
   use cint_tab_stg_w, only: STG_W
   implicit none
   public :: STG_X, STG_W
end module cint_tab_stg_roots
