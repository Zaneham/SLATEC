MODULE nonlin_eq
  USE service, ONLY : SP, DP
  USE lapack, ONLY : SHSEQR, CHSEQR  ! Required by rpqr79, cpqr79
  IMPLICIT NONE

CONTAINS
  include"cpevl.inc"
  include"cpevlr.inc"
  include"cpqr79.inc"
  include"cpzero.inc"
  include"dfzero.inc"
  include"dsos.inc"
  include"dsoseq.inc"
  include"dsossl.inc"
  include"fzero.inc"
  include"rpqr79.inc"
  include"rpzero.inc"
  include"sos.inc"
  include"soseqs.inc"
  include"sossol.inc"
END MODULE nonlin_eq