!> ZH: External linkage wrappers for LAPACK stubs
!> These provide the sgtsv_, sgetrs_, etc. symbols that gfortran expects
!> when implicit-external = true is set

SUBROUTINE sgtsv(N, Nrhs, Dl, D, Du, B, Ldb, Info)
  USE service, ONLY : SP
  INTEGER, INTENT(IN) :: Ldb, N, Nrhs
  INTEGER, INTENT(OUT) :: Info
  REAL(SP), INTENT(INOUT) :: B(Ldb,*), D(*), Dl(*), Du(*)
  Info = -999
END SUBROUTINE

SUBROUTINE dgtsv(N, Nrhs, Dl, D, Du, B, Ldb, Info)
  USE service, ONLY : DP
  INTEGER, INTENT(IN) :: Ldb, N, Nrhs
  INTEGER, INTENT(OUT) :: Info
  REAL(DP), INTENT(INOUT) :: B(Ldb,*), D(*), Dl(*), Du(*)
  Info = -999
END SUBROUTINE

SUBROUTINE sgetrs(Trans, N, Nrhs, A, Lda, Ipiv, B, Ldb, Info)
  USE service, ONLY : SP
  CHARACTER, INTENT(IN) :: Trans
  INTEGER, INTENT(IN) :: Lda, Ldb, N, Nrhs
  INTEGER, INTENT(OUT) :: Info
  INTEGER, INTENT(IN) :: Ipiv(*)
  REAL(SP), INTENT(IN) :: A(Lda,*)
  REAL(SP), INTENT(INOUT) :: B(Ldb,*)
  Info = -999
END SUBROUTINE

SUBROUTINE dgetrs(Trans, N, Nrhs, A, Lda, Ipiv, B, Ldb, Info)
  USE service, ONLY : DP
  CHARACTER, INTENT(IN) :: Trans
  INTEGER, INTENT(IN) :: Lda, Ldb, N, Nrhs
  INTEGER, INTENT(OUT) :: Info
  INTEGER, INTENT(IN) :: Ipiv(*)
  REAL(DP), INTENT(IN) :: A(Lda,*)
  REAL(DP), INTENT(INOUT) :: B(Ldb,*)
  Info = -999
END SUBROUTINE

SUBROUTINE cgetrs(Trans, N, Nrhs, A, Lda, Ipiv, B, Ldb, Info)
  USE service, ONLY : SP
  CHARACTER, INTENT(IN) :: Trans
  INTEGER, INTENT(IN) :: Lda, Ldb, N, Nrhs
  INTEGER, INTENT(OUT) :: Info
  INTEGER, INTENT(IN) :: Ipiv(*)
  COMPLEX(SP), INTENT(IN) :: A(Lda,*)
  COMPLEX(SP), INTENT(INOUT) :: B(Ldb,*)
  Info = -999
END SUBROUTINE

SUBROUTINE sgbtrs(Trans, N, Kl, Ku, Nrhs, Ab, Ldab, Ipiv, B, Ldb, Info)
  USE service, ONLY : SP
  CHARACTER, INTENT(IN) :: Trans
  INTEGER, INTENT(IN) :: Kl, Ku, Ldab, Ldb, N, Nrhs
  INTEGER, INTENT(OUT) :: Info
  INTEGER, INTENT(IN) :: Ipiv(*)
  REAL(SP), INTENT(IN) :: Ab(Ldab,*)
  REAL(SP), INTENT(INOUT) :: B(Ldb,*)
  Info = -999
END SUBROUTINE

SUBROUTINE dgbtrs(Trans, N, Kl, Ku, Nrhs, Ab, Ldab, Ipiv, B, Ldb, Info)
  USE service, ONLY : DP
  CHARACTER, INTENT(IN) :: Trans
  INTEGER, INTENT(IN) :: Kl, Ku, Ldab, Ldb, N, Nrhs
  INTEGER, INTENT(OUT) :: Info
  INTEGER, INTENT(IN) :: Ipiv(*)
  REAL(DP), INTENT(IN) :: Ab(Ldab,*)
  REAL(DP), INTENT(INOUT) :: B(Ldb,*)
  Info = -999
END SUBROUTINE

SUBROUTINE cgbtrs(Trans, N, Kl, Ku, Nrhs, Ab, Ldab, Ipiv, B, Ldb, Info)
  USE service, ONLY : SP
  CHARACTER, INTENT(IN) :: Trans
  INTEGER, INTENT(IN) :: Kl, Ku, Ldab, Ldb, N, Nrhs
  INTEGER, INTENT(OUT) :: Info
  INTEGER, INTENT(IN) :: Ipiv(*)
  COMPLEX(SP), INTENT(IN) :: Ab(Ldab,*)
  COMPLEX(SP), INTENT(INOUT) :: B(Ldb,*)
  Info = -999
END SUBROUTINE
