Subroutine finiteStrainNew(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use loadingLib
	use finiteStrainLib
	
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
	Integer :: NparamElem, constLaw 
	integer , parameter ::  NdimE = 2 , NodG = 3 ! this code only makes sense in 2D and for triangles 
	Real*8 :: PhiA,PhiB, dv , Det , bm(NdimE), MatPar(9)
	Real*8 :: Psi(NdimE,NGPMAX), W(NGPMAX), Phi(NodG), dPhiA(NdimE), dPhiB(NdimE),&
			  dPhi_G(NdimE,NodG), dPhi_L(NdimE,NodG), Jac(NDimE,NDimE)
	Real*8 :: GradU(NdimE,NdimE) , F(NdimE,NdimE)
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE) 
	Real*8 :: Sol1Mod(NodG*NdimE) , XLLMod(NodG*NdimE), DGradU(NdimE,NdimE)
	integer :: pOrder , NGP
	Real*8 :: LoadPar(6) , energy
	integer :: iShiftFluc , LoadType , LoadTypeProg


!~ 	write(*,*) "TLcompVolPeriodicLag2D"
	iShiftFluc = nint(CommonPar(1))
	constLaw = nint(CommonPar(2))
	matPar(1:3) = commonPar(3:5)
	pOrder = nint(commonPar(6))
	LoadType = nint(commonPar(7))
	LoadTypeProg = nint(commonPar(8))
	LoadPar(:) = CommonPar(9:14)
	
	call getSlice(Sol1Mod,Sol1,1,NodG,iShiftFluc + 1 ,iShiftFluc + NdimE, iDofT)
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
		
		call calcGradU(GradU,Sol1Mod,dPhi_G,NdimE,NodG,NdimE)
		call setF(F,LoadPar,Time,DelT,LoadTypeProg,LoadType)
		
		F = F + deltaKron(1:NdimE,1:NdimE) + GradU
		
		call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
		call calcD(D,F,NdimE,matPar,constLaw)	
		
		energy = 0.0d0
		
		DV=Det*W(nG)
				
		call T4xT2(DGradU,D,GradU)

		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftFluc
			PhiA = Phi(A)
			dPhiA=dPhi_G(:,A)
			
			Do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = BE(App) - dot_product(Ppk(p,:),dPhiA)*DV  				
				BE(App) = BE(App) + dot_product(DGradU(p,:),dphiA)*DV  
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*iDofT + iShiftFluc
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
 	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine finiteStrainNewS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftFluc 
    integer , parameter :: NodG = 3 , NdimE = 2

    Coupling = 0
    
    iShiftFluc = nint(commonPar(1))
    
    do A = 1,NodG
        ApRow = (A-1) * iDofT + iShiftFluc
		
		do B = 1,NodG
			BpCol = (B-1) * iDofT + iShiftFluc

			do p=1,nDimE
			do q=1,nDimE

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo
 	
	call numPrint(Coupling)

end Subroutine
