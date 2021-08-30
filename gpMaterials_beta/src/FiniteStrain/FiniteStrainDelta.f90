Subroutine finiteStrainDelta(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use loadingLib
	use finiteStrainLib
	use damageLib
	
	implicit none

	!	PARAMETERS
	integer :: iSimplex, NGPMAX, iBu 		
	parameter (iBu = 0, NGPMAX = 11, iSimplex = 0) ! No Bubbles and Linear, ! iSimplex = 0 -> The element is a Simplex 
!~ 	parameter (iBu = 0, pOrder = 2,NGP = 11, iSimplex = 0) ! No Bubbles and Linear, ! iSimplex = 0 -> The element is a Simplex 

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: NparamElem, constLaw ! if NparamDamage = 0 , no damage
	logical :: ifAnalytical
	integer , parameter ::  NdimE = 2 , NodG = 3 ! this code only makes sense in 2D and for triangles 
	Real*8 :: PhiA,PhiB, dv , Det , bm(NdimE), MatPar(9)
	Real*8 :: Psi(NdimE,NGPMAX), W(NGPMAX), Phi(NodG), dPhiA(NdimE), dPhiB(NdimE),&
			  dPhi_G(NdimE,NodG), dPhi_L(NdimE,NodG), Jac(NDimE,NDimE)
	Real*8 :: GradU(NdimE,NdimE) , GradDeltaU(NdimE,NdimE), F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), DB(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE)
	Real*8 :: SolU(NodG*NdimE) , SolDeltaU(NodG*NdimE), XLLMod(NodG*NdimE), DBGradDeltaU(NdimE,NdimE)
	integer :: pOrder , NGP
	Real*8 :: LoadPar(6) , energy
	integer :: iShiftU, iShiftDeltaU , LoadType , LoadTypeProg
	!! commonPar for the damage starts after the commonPar of the actual function
	integer , parameter :: iShiftDamageParam = 0 , ishiftDamageCommonPar = 16
	integer :: IdamageModify

	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	constLaw = nint(CommonPar(3))
	matPar(1:3) = commonPar(4:6)
	pOrder = nint(commonPar(7))
	LoadType = nint(commonPar(8))
	LoadTypeProg = nint(commonPar(9))
	LoadPar(:) = CommonPar(10:15)
	IdamageModify = nint(commonPar(16))
	
	call getSlice(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSlice(SolDeltaU,Sol1,1,NodG,iShiftDeltaU + 1 ,iShiftDeltaU + NdimE, iDofT)
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
			
	Do nG = 1, NGP ! LoopGauss
	
		Call LocalShapeDer(NodG,NdimE,dPhi_L,PSI(1,nG),pOrder,iBu) ! PSI(1,nG) is a pointer to the nG^th column
		if ( nG .eq. 1 ) then 
			Call Jacobian(Jac,Det,NdimE,NodG,XLLmod,dPhi_L)
		endif
		Call GlobalShapeDer(NodG,NdimE,dPhi_L,dPhi_G,Jac)
		Call ShapeF (NodG,NdimE,Phi,PSI(1,nG),pOrder,iBu)
		
		call calcGradU(GradU,SolU,dPhi_G,NdimE,NodG,NdimE)
		call calcGradU(GradDeltaU,SolDeltaU,dPhi_G,NdimE,NodG,NdimE)
		
		call setF(F,LoadPar,Time,DelT,LoadTypeProg,LoadType)
		
		F = F + deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		energy = 0.0d0
		if(IdamageModify == 1 .or. IdamageModify == 5 ) then ! Coupled
			call strainEnergy(energy,F,matPar,constLaw)
			call damageUpdate(commonPar,Param,energy, iShiftDamageCommonPar,iShiftDamageParam,nG)
		end if
		
		call damageModifyTangentNew(D,Ppk,commonPar,Param, NdimE, &
									iShiftDamageCommonPar,iShiftDamageParam,nG,IdamageModify)

		
!~ 		DB = 0.0d0
!~ 		
!~ 		call damageModifyTangentDouble(D,DB,Ppk,commonPar,Param,NdimE, &
!~ 									iShiftDamageCommonPar,iShiftDamageParam,IdamageModify)
!~ 		
		DV=Det*W(nG)
		
			
!~ 		call T4xT2(DBGradDeltaU,DB,GradDeltaU)

		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftDeltaU
			PhiA = Phi(A)
			dPhiA=dPhi_G(:,A)
			
			Do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = BE(App) - dot_product(Ppk(p,:),dPhiA)*DV  
!~ 				BE(App) = BE(App) + dot_product(DBGradDeltaU(p,:),dphiA)*DV  				  
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*iDofT + iShiftDeltaU
					PhiB   = Phi(B)
					dPhiB=dPhi_G(:,B) 
 					
					do q=1,NdimE
						Bqp = BpCol+q
						
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(D,dPhiB,dPhiA,p,q,NdimE)*DV  
				
					end do ! loop Bq
				Enddo !LoopCol
			end do ! loop Ap
		Enddo !LoopRow
	EndDo !LoopGauss
 	
 	
 	Do A=1, NodG !LoopRow
		ApRow  = (A-1)*iDofT + iShiftU
		
		Do p = 1 , NdimE
			App = ApRow+p
			
			BE(App) = SolU((A-1)*NdimE + p)  				  
			
			AE(App,App) = 1.0d0
			AE(App,App + iShiftDeltaU - iShiftU) = -1.0d0
			
		end do ! loop Ap
	Enddo !LoopRow
 	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine finiteStrainDeltaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftU , iShiftDeltaU
    integer , parameter :: NodG = 3 , NdimE = 2

    Coupling = 0
    
    iShiftU = nint(commonPar(1))
    iShiftDeltaU = nint(commonPar(2))
    
    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftDeltaU
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftDeltaU

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo


    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftU
		do p=1,nDimE
			Coupling (ApRow + p, ApRow + p) = 1
			Coupling (ApRow + p, ApRow + p + iShiftDeltaU - iShiftU) = 1
		enddo
	enddo
 	
	call numPrint(Coupling)

end Subroutine
