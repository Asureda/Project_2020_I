! Modul per les condicions periodiques de contorn
MODULE PBC
use READ_DATA
implicit none
contains

FUNCTION PBC1(x,L)
! Funcio que retorna la distancia donades unes condicions periodiques de contorn i la longitud de la caixa de simulacio L
    IMPLICIT NONE
    REAL*8 x,L,PBC1
    PBC1=x-int(2d0*x/L)*L
    RETURN
END FUNCTION

FUNCTION PBC2(x,L)
! Funcio que retorna la distancia donades unes condicions periodiques de contorn i la longitud de la caixa de simulacio L
    IMPLICIT NONE
    REAL*8 x,L,PBC2
    IF(x.lt.0)THEN
        x=L+mod(x,L)
    ELSE IF(x.gt.L) THEN
        x=mod(x,L)
    END IF
    PBC2=x
    RETURN
END FUNCTION

END MODULE PBC