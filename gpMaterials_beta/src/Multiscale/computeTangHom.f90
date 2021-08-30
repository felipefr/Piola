!     ------------------------------------------------------------------
Subroutine computeTangHom(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use loadingLib
	use finiteStrainLib !! New Way	
	use damageLib
	
	implicit none

	!	PARAMETERS
	integer :: iSimplex, NGPMAX, iBu 		
	parameter (iBu = 0, NGPMAX = 11, iSimplex = 0) ! No Bubbles and Linear, ! iSimplex = 0 -> The element is a Simplex 
!~ 	parameter (iBu = 0, pOrder = 2,NGP = 11, iSimplex = 0) ! No Bubbles and Linear, ! iSimplex = 0 -> The element is a Simplex 

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l, A, ApRow! counters 
	integer ::  ipDim, ipDimE, ipDofU,ipDofUmod, ipDofKL,ipDofKLmod , p, q
	Integer :: NparamElem, constLaw ! if NparamDamage = 0 , no damage
	logical :: ifAnalytical
	integer, parameter :: NdimE = 2, NodG = 3 ! this code only makes sense in 2D 
	Real*8 :: PhiA,PhiB, dv , Det , MatPar(9)
	Real*8 :: Psi(NdimE,NGPMAX), W(NGPMAX), Phi(NodG), dPhiA(NdimE), dPhiB(NdimE),&
			  dPhi_G(NdimE,NodG), dPhi_L(NdimE,NodG), Jac(NDimE,NDimE)
	Real*8 :: U(NdimE), GradU(NdimE,NdimE) , F(NdimE,NdimE), FT(NdimE,NdimE) , FinvT(NdimE,NdimE) , C(NdimE,NdimE) , detF
	Real*8 :: Dthom(NdimE,NdimE,NdimE,NdimE) , DtGradDeltaUkl(NdimE,NdimE)
	Real*8 :: Spk(NdimE,NdimE), Dt(NdimE,NdimE,NdimE,NdimE), Dc(NdimE,NdimE,NdimE,NdimE) , Ppk(NdimE,NdimE)
	Real*8 :: SolU(NodG*NdimE) , SolDeltaUkl(NodG*NdimE,NdimE,NdimE), XLLMod(NodG*NdimE), GradDeltaUkl(NdimE,NdimE)
	integer :: pOrder , NGP, IdamageModify
	integer :: iShiftDeltaU , iShiftU , iShiftC, iAux, ipDofDeltaU, NElem
	!! commonPar for the damage starts after the commonPar of the actual function
	integer , parameter :: iShiftDamageParam = 0 , ishiftDamageCommonPar = 10
	
	iShiftC = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iShiftU = nint(CommonPar(3))
	constLaw = nint(CommonPar(4))
	matPar(1:3) = commonPar(5:7)
	pOrder = nint(commonPar(8))
	NElem = nint(commonPar(9))
	IdamageModify = nint(commonPar(10)) 
	
	call getSlice(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	do k = 1, NdimE
	do l = 1, NdimE
		call getSlice(SolDeltaUkl(:,k,l),Sol1,1,NodG,iShiftDeltaU + iShiftKL(k,l,NdimE)*NdimE + 1 , &
		             iShiftDeltaU + iShiftKL(k,l,NdimE)*NdimE +  NdimE, iDofT)
	end do
	end do
	call getSlice(XLLmod,XLL,1,NodG,1 ,NdimE, Ndim)

!~ 	
	select case(pOrder)
		case(1)
			NGP = 1
		case(2)
			NGP = 3
		case default
			write(0,*) "pOrder must be 1 or 2"
	end select
	
	W = 0.0d0		
	Call GaussRule (Psi,W,NdimE,NGP,iSimplex)
	
	AE = 0.0d0
	BE = 0.0d0
	
	!--------------------------
	Do nG = 1, NGP ! LoopGauss
	! dPhi_L = Phirst derivatives with respect to intrinsic Local coordinates refers to linear shape functions
	! dPhi_G = Phirst derivatives with respect to Global coordinates
		Call LocalShapeDer(NodG,NdimE,dPhi_L,PSI(1,nG),pOrder,iBu) ! PSI(1,nG) is a pointer to the nG^th column
		if ( nG .eq. 1 ) then 
			Call Jacobian(Jac,Det,NdimE,NodG,XLLmod,dPhi_L)
		endif
		Call GlobalShapeDer(NodG,NdimE,dPhi_L,dPhi_G,Jac)
		Call ShapeF (NodG,NdimE,Phi,PSI(1,nG),pOrder,iBu)
				
		!!! New way
		call calcGradU(GradU,SolU,dPhi_G,NdimE,NodG,NdimE)
		F = deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(Dt,F,NdimE,matPar,constLaw)	
		
		call damageModifyTangentNew(Dt,Ppk,commonPar,Param,NdimE, &
									iShiftDamageCommonPar,iShiftDamageParam, nG, IdamageModify)
		
		DV=Det*W(nG)
						
		Dthom = Dt
 		
		do k = 1, NdimE
		do l = 1, NdimE
			call calcGradU(GradDeltaUkl,SolDeltaUkl(:,k,l),dPhi_G,NdimE,NodG,NdimE)
			call T4xT2(DtGradDeltaUKl,Dt,GradDeltaUkl)
			Dthom(:,:,k,l) = Dthom(:,:,k,l) + DtGradDeltaUkl  
		end do
		end do 
		
		A = NodG + 1
		ApRow  = (A-1)*iDofT + iShiftC
			
		do i = 1 , NdimE
		do j = 1 , NdimE
		do k = 1 , NdimE
		do l = 1 , NdimE
			iAux = ApRow + iShiftIJKL(i,j,k,l,NdimE) + 1
			BE(iAux) = Dthom(i,j,k,l)*dV
			AE(iAux,iAux) = 1.0d0/real(NElem) 
		end do
		end do
		end do
		end do
!~ 	
		BE(ApRow + 5) = Ppk(1,1)*dV
		BE(ApRow + 6) = Ppk(1,2)*dV
		BE(ApRow + 7) = Ppk(2,1)*dV
		BE(ApRow + 8) = Ppk(2,2)*dV
	
	end do !LoopGauss

end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine computeTangHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	use funcAux
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A, ApRow, iShiftC,  i, j, k , l, iAux
	integer, parameter :: NodG = 3, NdimE = 2

    Coupling = 0
    
    iShiftC = nint(commonPar(1))

	A = NodG + 1
	ApRow  = (A-1)*iDofT + iShiftC
		
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		iAux = ApRow + iShiftIJKL(i,j,k,l,NdimE) + 1
		Coupling(iAux,iAux) = 1
	end do
	end do
	end do
	end do
		
!~ 	write(*,*) "TLflucTangHom Symbolic"
!~ 	call IprintMat(Coupling)
end Subroutine
