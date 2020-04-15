!GRUP I: Àlex, Oriol, Laia, Sílvia i Elena
module PBC

use READ_DATA

implicit none

contains

FUNCTION PBC1(x,L)

    !OBJECTIU: Definir les condicions periòdiques de contorn per tal de tenir les partícules
    !           definides sempre en una caixa i no variar així el nombre de partícules de la nostra simulació.
    
    !INPUTS: posicions(x), costat de la caixa (L)
    
    !OUTPUTS: posicions(x)
    
    IMPLICIT NONE
    
    REAL*8 x,L,PBC1
    
    PBC1=x-int(2d0*x/L)*L
    
    RETURN
END FUNCTION

FUNCTION PBC2(x,L)

    !OBJECTIU: Definir les condicions periòdiques de contorn per tal de tenir les partícules
    !           definides sempre en una caixa i no variar així el nombre de partícules de la nostra simulació.
    
    !INPUTS: posicions(x), costat de la caixa (L)
    
    !OUTPUTS: posicions(x)
    
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
