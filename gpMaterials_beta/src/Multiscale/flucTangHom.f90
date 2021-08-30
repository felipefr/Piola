!     ------------------------------------------------------------------
Subroutine flucTangHom(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use finiteStrainLib
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
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw ! if NparamDamage = 0 , no damage
	Integer , parameter :: NdimE = 2, NodG = 3 ! this code only makes sense in 2D 
	Real*8 :: PhiA,PhiB, dv , Det , MatPar(9)
	Real*8 :: Psi(NdimE,NGPMAX), W(NGPMAX), Phi(NodG), dPhiA(NdimE), dPhiB(NdimE),&
			  dPhi_G(NdimE,NodG), dPhi_L(NdimE,NodG), Jac(NDimE,NDimE)
	Real*8 :: U(NdimE), GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: DGradU(NdimE,NdimE), DFcan(NdimE,NdimE)
	Real*8 :: Ppk(NdimE,NdimE), D(NdimE,NdimE,NdimE,NdimE)
	Real*8 :: Sol1Mod(NodG*NdimE) , XLLMod(NodG*NdimE), Fcan(NdimE,NdimE) 
!~ 	Real*8 :: DDGradPhiBdotGradPhiA , dot_product2 ! function
	integer :: pOrder , NGP
	integer :: iShiftFlucKL , iShiftU, kl, IdamageModify
	!! commonPar for the damage starts after the commonPar of the actual function
	integer , parameter :: iShiftDamageParam = 0 , ishiftDamageCommonPar = 9
	
!~ 	write(0,*) "flucTangHom"
	
	iShiftFlucKL = nint(CommonPar(1))
	iShiftU = nint(CommonPar(2))
	kl = nint(CommonPar(3))
	constLaw = nint(CommonPar(4))
	matPar(1:3) = commonPar(5:7)
	pOrder = nint(commonPar(8))
	IdamageModify = nint(commonPar(9))  
	
	call getIJfromK(k,l,kl)
	
	Fcan = 0.0d0
	Fcan(k,l) = 1.0d0
			
	call getSlice(Sol1Mod,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSlice(XLLmod,XLL,1,NodG,1 ,NdimE, Ndim)

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
		
	    call calcGradU(GradU,Sol1Mod,dPhi_G,NdimE,NodG,NdimE)
		
		F =  deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		call damageModifyTangentNew(D,Ppk,commonPar,Param, NdimE, &
								iShiftDamageCommonPar,iShiftDamageParam,nG,IdamageModify)
 		
		DV=Det*W(nG)
	
		DFcan = 0.0d0		
		call T4xT2(DFcan,D,Fcan)

		do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftFlucKL
			PhiA = Phi(A)
			dPhiA=dPhi_G(:,A)
			
			do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = BE(App) - dot_product(DFcan(p,:),dphiA)*DV  
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*iDofT + iShiftFlucKL
					PhiB   = Phi(B)
					dPhiB=dPhi_G(:,B) 
 					
					do q=1,NdimE
						Bqp = BpCol+q
						
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(D,dPhiB,dPhiA,p,q,NdimE)*DV  
				
					end do ! loop Bq
				end do !LoopCol
			end do ! loop Ap
		end do !LoopRow
	end do !LoopGauss
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine flucTangHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol, iShiftFlucKL
    integer , parameter :: NodG = 3, NdimE = 2 
   
    Coupling = 0
    
    iShiftFlucKL = nint(commonPar(1))
    
    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftFlucKL
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftFlucKL

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo
	
!~ 	write(*,*) "TLflucTangHom Symbolic"
!~ 	call IprintMat(Coupling)

end Subroutine
