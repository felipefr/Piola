Subroutine FbarSpatial(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use loadingLib
	use finiteStrainLib
	use ptsGaussLib2
	
	implicit none
!~ 
!~ 	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ipDim, ipDimE, ipDof,ipDofmod ! counters 
	Integer :: constLaw 
	integer , parameter ::  NdimE = 2! this code only makes sense in 2D and for triangles 
	Real*8 :: PhiA,PhiB, dV , Det , bm(NdimE), MatPar(9)
	Real*8 :: dPhiA(NdimE), dPhiB(NdimE), dPhiB0(NdimE), detFs0, detFs
	Real*8 , dimension(NdimE,NdimE) :: GradUs , GradUs0, DGradUs, sigma, Fs, Fs0, Fbar
	Real*8 , dimension(NdimE,NdimE,NdimE,NdimE) :: Ds , tenQs
	real*8 , allocatable ::  SolU(:) , Xel(:) !! all have dimension NodG*NdimE
	integer :: pOrder , NGP, NodG, iFEMtype , iSimplex, iBubble
	Real*8 :: LoadPar(6) , he, energy
	integer :: iShiftU, iShiftDeltaU , LoadType , LoadTypeProg
	type(ptGaussClass) :: PtG, PtG0
	
	
	iShiftU = nint(CommonPar(1))
	iShiftDeltaU = nint(CommonPar(2))
	iFEMtype =  nint(commonPar(3)) 
	constLaw = nint(CommonPar(4))
	matPar(1:3) = commonPar(5:7)
	LoadType = nint(commonPar(8))
	LoadTypeProg = nint(commonPar(9))
	LoadPar(:) = CommonPar(10:15)
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)

	call getSliceAllocate(SolU,Sol1,1,NodG,iShiftU + 1 ,iShiftU + NdimE, iDofT)
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
	
	Xel = Xel + SolU
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
		
	nG = 1
	call PtG0%init(Xel,NodG,NdimE,1,pOrder,iBubble,iSimplex)
	call PtG0%calcGradU(GradUs0,SolU,nG)
	
	call setF(Fs0,LoadPar,Time,DelT,LoadTypeProg,LoadType)

	call computeFs(Fs0,detFs0,gradUs0,NdimE)
			
	Do nG = 1, NGP ! LoopGauss
	
		GradUs = 0.0d0 
		
		call PtG%calcGradU(GradUs,SolU,nG)
		
		call setF(Fs,LoadPar,Time,DelT,LoadTypeProg,LoadType)
		
		call computeFs(Fs,detFs,gradUs,NdimE)

		Fbar = ((detFs0/detFs)**(1.0d0/2.0d0))*Fs !! for 2D
		
		call calcPpk(sigma,Fbar,NdimE,matPar,constLaw)		
		call calcD(Ds,Fbar,NdimE,matPar,constLaw)	
		call modifyToSpatial(Ds,sigma,Fbar,detFs0,NdimE) 
			
		call build_tenQs(Ds,sigma,tenQs,NdimE)
		
		dV=ptG%dV(nG)
		
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftDeltaU
			PhiA = PtG%Phi(A,nG)
			dPhiA= PtG%dPhi_G(:,A,nG)
			
			Do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = BE(App) - dot_product(sigma(p,:),dPhiA)*DV  
				
				do B= 1 , NodG ! LoopCol ! not considering the simmetry
					BpCol  = (B-1)*iDofT + iShiftDeltaU
					PhiB  = PtG%Phi(B,nG)
					dPhiB = PtG%dPhi_G(:,B,nG) 
					dPhiB0 = PtG0%dPhi_G(:,B,1)
			
					do q=1,NdimE
						Bqp = BpCol+q
						
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(Ds,dPhiB,dPhiA,p,q,NdimE)*DV 
						AE(App,Bqp) = AE(App,Bqp) + DGradPhiBdotGradPhiA(tenQs,dPhiB0 - dPhiB,dPhiA,p,q,NdimE)*DV
!~ 						write(0,*) AE(App,Bqp)
					end do ! loop Bq
				Enddo !LoopCol
			end do ! loop Ap
		Enddo !LoopRow
	EndDo !LoopGauss
 	
 	if(iShiftU /= iShiftDeltaU) then
		Do A=1, NodG !LoopRow
			ApRow  = (A-1)*iDofT + iShiftU
			
			Do p = 1 , NdimE
				App = ApRow+p
				
				BE(App) = SolU((A-1)*NdimE + p)  				  
				
				AE(App,App) = 1.0d0
				AE(App,App + iShiftDeltaU - iShiftU) = -1.0d0
				
			end do ! loop Ap
		end do !LoopRow
	end if
	
end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine FbarSpatialS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    use ptsGaussLib2
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol,iShiftU , iShiftDeltaU, Femtype, NodG
    integer , parameter :: NdimE = 2

    Coupling = 0
    
    iShiftU = nint(commonPar(1))
    iShiftDeltaU = nint(commonPar(2))
    femType = nint(commonPar(3)) 
    
    call setNodG(femtype, NodG)

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

 	if(iShiftU /= iShiftDeltaU) then
		do A = 1,NodG
			ApRow = (A-1) * iDofT + iShiftU
			do p=1,nDimE
				Coupling (ApRow + p, ApRow + p) = 1
				Coupling (ApRow + p, ApRow + p + iShiftDeltaU - iShiftU) = 1
			enddo
		enddo
	end if

 	
	call numPrint(Coupling)

end Subroutine

subroutine computeFs(Fs,detFs,gradUs,NdimE)
	use funcAux
	implicit none
	
	real*8 , intent(inout) :: Fs(NdimE,NdimE) 
	real*8, intent(out) :: detFs
	real*8, intent(in) :: gradUs(NdimE,NdimE)
	integer , intent(in) :: NdimE
	
	real*8 :: Fsinv(NdimE,NdimE), detFdummy
	
	Fs = Fs + deltaKron(1:NdimE,1:NdimE) !! Fs can be prescribed
	
	call matInv(FsInv,detFdummy,Fs)
	FsInv = Fsinv - gradUs
	
	call matInv(Fs,detFs,FsInv)
	detFs = 1.0d0/detFs

end subroutine




