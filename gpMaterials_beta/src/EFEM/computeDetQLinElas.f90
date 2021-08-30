Subroutine computeMinDetQLinElas &
    (AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

	use funcAux
	use multiscaleLib
	use damageNewLib
	
    implicit none
    
    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======
    integer, parameter :: NdimE = 2, nStepTheta = 999
    integer :: ip1 , ip2, iShiftC, iShiftMinDetQ, iShiftThetaCrit
    Real*8 ::  MinDetQ, ThetaCrit, betaCrit

	iShiftMinDetQ = nint(commonPar(1))
	iShiftThetaCrit = nint(commonPar(2))
	

!~ 	call getMinDetQ(minDetQ,ThetaCrit,betaCrit, Sol1,NdimE,idofT,1,iShiftC,nStepTheta) !! is '1' because is just one node
	
	
	Param(LengthParam + iShiftMinDetQ + 1) = 1.0d0
	Param(LengthParam + iShiftThetaCrit + 1) = 2.0d0
	
	AE(1,1) = 1.0d0
	BE(1) =  1.0d0
		
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeMinDetQLinElasS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
	
    implicit none
    
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer ::  Coupling(MaxLRows,MaxLRows)
    Real*8 :: CommonPar(*)
    integer :: iShiftMinDetQ, iShiftThetaCrit, ip1, ip2
	
	iShiftMinDetQ = nint(commonPar(1))
	iShiftThetaCrit = nint(commonPar(2))

	Coupling(1,1) = 1
	
end Subroutine

