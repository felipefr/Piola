MODULE DETERMINANT
! Used in matrixcoef and in some elements and in writeElemOut, send to elements the specified variables
! Verified if used in jacobian
        REAL*8  :: DET
        INTEGER :: LengthParam,LengthJParam,iElem,IdElemLib
        INTEGER :: iWriteElemUnit,Id_Write, OUnit_residualSensitivity , OUnit_residualStress
END MODULE DETERMINANT
