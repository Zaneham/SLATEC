MODULE data_handling
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  include"c1merg.inc"
  include"d1merg.inc"
  include"dpperm.inc"
  include"dpsort.inc"
  include"dsort.inc"
  include"hpperm.inc"
  include"hpsort.inc"
  include"i1merg.inc"
  include"ipperm.inc"
  include"ipsort.inc"
  include"isort.inc"
  include"s1merg.inc"
  include"spperm.inc"
  include"spsort.inc"
  include"ssort.inc"
END MODULE data_handling