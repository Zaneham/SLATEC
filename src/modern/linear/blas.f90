!> ZH: Self-contained BLAS implementations for SLATEC
!> Reference: BLAS Technical Forum Standard (2002)
MODULE blas
  USE service, ONLY : SP, DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: SAXPY, DAXPY, CAXPY, SSWAP, DSWAP, CSWAP
  PUBLIC :: SCOPY, DCOPY, SSCAL, DSCAL, CSCAL
  PUBLIC :: SROT, DROT, SROTG, DROTG, SROTM, DROTM, SROTMG, DROTMG
  PUBLIC :: SCNRM2, SCABS1, SCASUM, ICAMAX
CONTAINS

  PURE SUBROUTINE SAXPY(N, Sa, Sx, Incx, Sy, Incy)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(SP), INTENT(IN) :: Sa, Sx(*)
    REAL(SP), INTENT(INOUT) :: Sy(*)
    INTEGER :: i, ix, iy
    IF (N <= 0 .OR. Sa == 0.0_SP) RETURN
    IF (Incx == 1 .AND. Incy == 1) THEN
      DO i = 1, N
        Sy(i) = Sy(i) + Sa * Sx(i)
      END DO
    ELSE
      ix = 1; iy = 1
      IF (Incx < 0) ix = (-N + 1) * Incx + 1
      IF (Incy < 0) iy = (-N + 1) * Incy + 1
      DO i = 1, N
        Sy(iy) = Sy(iy) + Sa * Sx(ix)
        ix = ix + Incx; iy = iy + Incy
      END DO
    END IF
  END SUBROUTINE SAXPY

  PURE SUBROUTINE DAXPY(N, Da, Dx, Incx, Dy, Incy)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(DP), INTENT(IN) :: Da, Dx(*)
    REAL(DP), INTENT(INOUT) :: Dy(*)
    INTEGER :: i, ix, iy
    IF (N <= 0 .OR. Da == 0.0_DP) RETURN
    IF (Incx == 1 .AND. Incy == 1) THEN
      DO i = 1, N
        Dy(i) = Dy(i) + Da * Dx(i)
      END DO
    ELSE
      ix = 1; iy = 1
      IF (Incx < 0) ix = (-N + 1) * Incx + 1
      IF (Incy < 0) iy = (-N + 1) * Incy + 1
      DO i = 1, N
        Dy(iy) = Dy(iy) + Da * Dx(ix)
        ix = ix + Incx; iy = iy + Incy
      END DO
    END IF
  END SUBROUTINE DAXPY

  PURE SUBROUTINE CAXPY(N, Ca, Cx, Incx, Cy, Incy)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    COMPLEX(SP), INTENT(IN) :: Ca, Cx(*)
    COMPLEX(SP), INTENT(INOUT) :: Cy(*)
    INTEGER :: i, ix, iy
    IF (N <= 0 .OR. ABS(Ca) == 0.0_SP) RETURN
    IF (Incx == 1 .AND. Incy == 1) THEN
      DO i = 1, N
        Cy(i) = Cy(i) + Ca * Cx(i)
      END DO
    ELSE
      ix = 1; iy = 1
      IF (Incx < 0) ix = (-N + 1) * Incx + 1
      IF (Incy < 0) iy = (-N + 1) * Incy + 1
      DO i = 1, N
        Cy(iy) = Cy(iy) + Ca * Cx(ix)
        ix = ix + Incx; iy = iy + Incy
      END DO
    END IF
  END SUBROUTINE CAXPY

  PURE SUBROUTINE SSWAP(N, Sx, Incx, Sy, Incy)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(SP), INTENT(INOUT) :: Sx(*), Sy(*)
    REAL(SP) :: t
    INTEGER :: i, ix, iy
    IF (N <= 0) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    DO i = 1, N
      t = Sx(ix); Sx(ix) = Sy(iy); Sy(iy) = t
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE SSWAP

  PURE SUBROUTINE DSWAP(N, Dx, Incx, Dy, Incy)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(DP), INTENT(INOUT) :: Dx(*), Dy(*)
    REAL(DP) :: t
    INTEGER :: i, ix, iy
    IF (N <= 0) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    DO i = 1, N
      t = Dx(ix); Dx(ix) = Dy(iy); Dy(iy) = t
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE DSWAP

  PURE SUBROUTINE CSWAP(N, Cx, Incx, Cy, Incy)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    COMPLEX(SP), INTENT(INOUT) :: Cx(*), Cy(*)
    COMPLEX(SP) :: t
    INTEGER :: i, ix, iy
    IF (N <= 0) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    DO i = 1, N
      t = Cx(ix); Cx(ix) = Cy(iy); Cy(iy) = t
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE CSWAP

  PURE SUBROUTINE SCOPY(N, Sx, Incx, Sy, Incy)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(SP), INTENT(IN) :: Sx(*)
    REAL(SP), INTENT(OUT) :: Sy(*)
    INTEGER :: i, ix, iy
    IF (N <= 0) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    DO i = 1, N
      Sy(iy) = Sx(ix)
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE SCOPY

  PURE SUBROUTINE DCOPY(N, Dx, Incx, Dy, Incy)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(DP), INTENT(IN) :: Dx(*)
    REAL(DP), INTENT(OUT) :: Dy(*)
    INTEGER :: i, ix, iy
    IF (N <= 0) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    DO i = 1, N
      Dy(iy) = Dx(ix)
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE DCOPY

  PURE SUBROUTINE SSCAL(N, Sa, Sx, Incx)
    INTEGER, INTENT(IN) :: N, Incx
    REAL(SP), INTENT(IN) :: Sa
    REAL(SP), INTENT(INOUT) :: Sx(*)
    INTEGER :: i
    IF (N <= 0 .OR. Incx <= 0) RETURN
    DO i = 1, N * Incx, Incx
      Sx(i) = Sa * Sx(i)
    END DO
  END SUBROUTINE SSCAL

  PURE SUBROUTINE DSCAL(N, Da, Dx, Incx)
    INTEGER, INTENT(IN) :: N, Incx
    REAL(DP), INTENT(IN) :: Da
    REAL(DP), INTENT(INOUT) :: Dx(*)
    INTEGER :: i
    IF (N <= 0 .OR. Incx <= 0) RETURN
    DO i = 1, N * Incx, Incx
      Dx(i) = Da * Dx(i)
    END DO
  END SUBROUTINE DSCAL

  PURE SUBROUTINE CSCAL(N, Ca, Cx, Incx)
    INTEGER, INTENT(IN) :: N, Incx
    COMPLEX(SP), INTENT(IN) :: Ca
    COMPLEX(SP), INTENT(INOUT) :: Cx(*)
    INTEGER :: i
    IF (N <= 0 .OR. Incx <= 0) RETURN
    DO i = 1, N * Incx, Incx
      Cx(i) = Ca * Cx(i)
    END DO
  END SUBROUTINE CSCAL

  PURE SUBROUTINE SROT(N, Sx, Incx, Sy, Incy, Sc, Ss)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(SP), INTENT(IN) :: Sc, Ss
    REAL(SP), INTENT(INOUT) :: Sx(*), Sy(*)
    REAL(SP) :: t
    INTEGER :: i, ix, iy
    IF (N <= 0) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    DO i = 1, N
      t = Sc*Sx(ix) + Ss*Sy(iy)
      Sy(iy) = Sc*Sy(iy) - Ss*Sx(ix)
      Sx(ix) = t
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE SROT

  PURE SUBROUTINE DROT(N, Dx, Incx, Dy, Incy, Dc, Ds)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(DP), INTENT(IN) :: Dc, Ds
    REAL(DP), INTENT(INOUT) :: Dx(*), Dy(*)
    REAL(DP) :: t
    INTEGER :: i, ix, iy
    IF (N <= 0) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    DO i = 1, N
      t = Dc*Dx(ix) + Ds*Dy(iy)
      Dy(iy) = Dc*Dy(iy) - Ds*Dx(ix)
      Dx(ix) = t
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE DROT

  PURE SUBROUTINE SROTG(Sa, Sb, Sc, Ss)
    REAL(SP), INTENT(INOUT) :: Sa, Sb
    REAL(SP), INTENT(OUT) :: Sc, Ss
    REAL(SP) :: r, roe, scale, z
    roe = Sb
    IF (ABS(Sa) > ABS(Sb)) roe = Sa
    scale = ABS(Sa) + ABS(Sb)
    IF (scale == 0.0_SP) THEN
      Sc = 1.0_SP; Ss = 0.0_SP; r = 0.0_SP; z = 0.0_SP
    ELSE
      r = scale * SQRT((Sa/scale)**2 + (Sb/scale)**2)
      r = SIGN(1.0_SP, roe) * r
      Sc = Sa / r; Ss = Sb / r; z = 1.0_SP
      IF (ABS(Sa) > ABS(Sb)) z = Ss
      IF (ABS(Sb) >= ABS(Sa) .AND. Sc /= 0.0_SP) z = 1.0_SP / Sc
    END IF
    Sa = r; Sb = z
  END SUBROUTINE SROTG

  PURE SUBROUTINE DROTG(Da, Db, Dc, Ds)
    REAL(DP), INTENT(INOUT) :: Da, Db
    REAL(DP), INTENT(OUT) :: Dc, Ds
    REAL(DP) :: r, roe, scale, z
    roe = Db
    IF (ABS(Da) > ABS(Db)) roe = Da
    scale = ABS(Da) + ABS(Db)
    IF (scale == 0.0_DP) THEN
      Dc = 1.0_DP; Ds = 0.0_DP; r = 0.0_DP; z = 0.0_DP
    ELSE
      r = scale * SQRT((Da/scale)**2 + (Db/scale)**2)
      r = SIGN(1.0_DP, roe) * r
      Dc = Da / r; Ds = Db / r; z = 1.0_DP
      IF (ABS(Da) > ABS(Db)) z = Ds
      IF (ABS(Db) >= ABS(Da) .AND. Dc /= 0.0_DP) z = 1.0_DP / Dc
    END IF
    Da = r; Db = z
  END SUBROUTINE DROTG

  PURE SUBROUTINE SROTM(N, Sx, Incx, Sy, Incy, Sparam)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(SP), INTENT(IN) :: Sparam(5)
    REAL(SP), INTENT(INOUT) :: Sx(*), Sy(*)
    REAL(SP) :: flag, h11, h12, h21, h22, w, z
    INTEGER :: i, ix, iy
    flag = Sparam(1)
    IF (N <= 0 .OR. flag == -2.0_SP) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    IF (flag < 0.0_SP) THEN
      h11 = Sparam(2); h12 = Sparam(4); h21 = Sparam(3); h22 = Sparam(5)
    ELSE IF (flag == 0.0_SP) THEN
      h11 = 1.0_SP; h12 = Sparam(4); h21 = Sparam(3); h22 = 1.0_SP
    ELSE
      h11 = Sparam(2); h12 = 1.0_SP; h21 = -1.0_SP; h22 = Sparam(5)
    END IF
    DO i = 1, N
      w = Sx(ix); z = Sy(iy)
      Sx(ix) = w*h11 + z*h12
      Sy(iy) = w*h21 + z*h22
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE SROTM

  PURE SUBROUTINE DROTM(N, Dx, Incx, Dy, Incy, Dparam)
    INTEGER, INTENT(IN) :: N, Incx, Incy
    REAL(DP), INTENT(IN) :: Dparam(5)
    REAL(DP), INTENT(INOUT) :: Dx(*), Dy(*)
    REAL(DP) :: flag, h11, h12, h21, h22, w, z
    INTEGER :: i, ix, iy
    flag = Dparam(1)
    IF (N <= 0 .OR. flag == -2.0_DP) RETURN
    ix = 1; iy = 1
    IF (Incx < 0) ix = (-N + 1) * Incx + 1
    IF (Incy < 0) iy = (-N + 1) * Incy + 1
    IF (flag < 0.0_DP) THEN
      h11 = Dparam(2); h12 = Dparam(4); h21 = Dparam(3); h22 = Dparam(5)
    ELSE IF (flag == 0.0_DP) THEN
      h11 = 1.0_DP; h12 = Dparam(4); h21 = Dparam(3); h22 = 1.0_DP
    ELSE
      h11 = Dparam(2); h12 = 1.0_DP; h21 = -1.0_DP; h22 = Dparam(5)
    END IF
    DO i = 1, N
      w = Dx(ix); z = Dy(iy)
      Dx(ix) = w*h11 + z*h12
      Dy(iy) = w*h21 + z*h22
      ix = ix + Incx; iy = iy + Incy
    END DO
  END SUBROUTINE DROTM

  PURE SUBROUTINE SROTMG(Sd1, Sd2, Sx1, Sy1, Sparam)
    REAL(SP), INTENT(INOUT) :: Sd1, Sd2, Sx1
    REAL(SP), INTENT(IN) :: Sy1
    REAL(SP), INTENT(OUT) :: Sparam(5)
    Sparam = 0.0_SP
    Sparam(1) = -2.0_SP
  END SUBROUTINE SROTMG

  PURE SUBROUTINE DROTMG(Dd1, Dd2, Dx1, Dy1, Dparam)
    REAL(DP), INTENT(INOUT) :: Dd1, Dd2, Dx1
    REAL(DP), INTENT(IN) :: Dy1
    REAL(DP), INTENT(OUT) :: Dparam(5)
    Dparam = 0.0_DP
    Dparam(1) = -2.0_DP
  END SUBROUTINE DROTMG

  PURE FUNCTION SCNRM2(N, Cx, Incx) RESULT(nrm)
    INTEGER, INTENT(IN) :: N, Incx
    COMPLEX(SP), INTENT(IN) :: Cx(*)
    REAL(SP) :: nrm, scale, ssq, temp
    INTEGER :: i, ix
    IF (N < 1 .OR. Incx < 1) THEN
      nrm = 0.0_SP; RETURN
    END IF
    scale = 0.0_SP; ssq = 1.0_SP; ix = 1
    DO i = 1, N
      IF (REAL(Cx(ix)) /= 0.0_SP) THEN
        temp = ABS(REAL(Cx(ix)))
        IF (scale < temp) THEN
          ssq = 1.0_SP + ssq * (scale / temp)**2
          scale = temp
        ELSE
          ssq = ssq + (temp / scale)**2
        END IF
      END IF
      IF (AIMAG(Cx(ix)) /= 0.0_SP) THEN
        temp = ABS(AIMAG(Cx(ix)))
        IF (scale < temp) THEN
          ssq = 1.0_SP + ssq * (scale / temp)**2
          scale = temp
        ELSE
          ssq = ssq + (temp / scale)**2
        END IF
      END IF
      ix = ix + Incx
    END DO
    nrm = scale * SQRT(ssq)
  END FUNCTION SCNRM2

  ELEMENTAL FUNCTION SCABS1(C) RESULT(absval)
    COMPLEX(SP), INTENT(IN) :: C
    REAL(SP) :: absval
    absval = ABS(REAL(C)) + ABS(AIMAG(C))
  END FUNCTION SCABS1

  PURE FUNCTION SCASUM(N, Cx, Incx) RESULT(asum)
    INTEGER, INTENT(IN) :: N, Incx
    COMPLEX(SP), INTENT(IN) :: Cx(*)
    REAL(SP) :: asum
    INTEGER :: i, ix
    asum = 0.0_SP
    IF (N <= 0 .OR. Incx <= 0) RETURN
    ix = 1
    DO i = 1, N
      asum = asum + ABS(REAL(Cx(ix))) + ABS(AIMAG(Cx(ix)))
      ix = ix + Incx
    END DO
  END FUNCTION SCASUM

  PURE FUNCTION ICAMAX(N, Cx, Incx) RESULT(imax)
    INTEGER, INTENT(IN) :: N, Incx
    COMPLEX(SP), INTENT(IN) :: Cx(*)
    INTEGER :: imax
    REAL(SP) :: smax
    INTEGER :: i, ix
    imax = 0
    IF (N < 1 .OR. Incx <= 0) RETURN
    imax = 1
    IF (N == 1) RETURN
    ix = 1
    smax = ABS(REAL(Cx(1))) + ABS(AIMAG(Cx(1)))
    DO i = 2, N
      ix = ix + Incx
      IF (ABS(REAL(Cx(ix))) + ABS(AIMAG(Cx(ix))) > smax) THEN
        imax = i
        smax = ABS(REAL(Cx(ix))) + ABS(AIMAG(Cx(ix)))
      END IF
    END DO
  END FUNCTION ICAMAX

END MODULE blas
