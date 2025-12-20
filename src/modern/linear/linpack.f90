include"dqdot/mpcom.inc"
include"slap/sslblk.inc"
include"slap/dslblk.inc"

MODULE linpack
  USE service, ONLY : SP, DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: CGBFA, CGEFA, CGECO, CPOCO, CPOSL, CGESL, CPOFA
  PUBLIC :: SGBFA, SGEFA, SGECO, SPOCO, SPOSL, SGESL, SGBSL, SPOFA
  PUBLIC :: DGBFA, DGEFA, DGECO, DPOCO, DPOSL, DGESL, DGBSL

CONTAINS
  include"linpack/csign1.inc"
  include"linpack/cgbfa.inc"
  include"linpack/cgeco.inc"
  include"linpack/cgefa.inc"
  include"linpack/cgesl.inc"
  include"linpack/cpoco.inc"
  include"linpack/cpofa.inc"
  include"linpack/cposl.inc"
  include"linpack/dgbfa.inc"
  include"linpack/dgbsl.inc"
  include"linpack/dgeco.inc"
  include"linpack/dgefa.inc"
  include"linpack/dgesl.inc"
  include"linpack/dpoco.inc"
  include"linpack/dpofa.inc"
  include"linpack/dposl.inc"
  include"linpack/sgbfa.inc"
  include"linpack/sgbsl.inc"
  include"linpack/sgeco.inc"
  include"linpack/sgefa.inc"
  include"linpack/sgesl.inc"
  include"linpack/spoco.inc"
  include"linpack/spofa.inc"
  include"linpack/sposl.inc"
END MODULE linpack
