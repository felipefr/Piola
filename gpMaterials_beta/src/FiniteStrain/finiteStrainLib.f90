module finiteStrainLib

use funcAux
use constitutiveLib

implicit none

public calcGradU , calcU, DGradPhiBdotGradPhiA, DGradGrad, calcPpk, calcD, modifyToSpatial, calcFC,&
		VirtualPower, buildFbarTensors, build_tenQs

contains

subroutine getPKmuFromMatrix(PKmu,Vmu,SolUf,Xel,matpar,Param,constLaw,iFtype,NdimE,iFemType)
	use damageNewLib
	use ptsGaussLib2
	use globalVariables , only : getF

	implicit none
		
	real*8 , intent(out) :: PKmu(:,:) , Vmu
	real*8, intent(in) :: matPar(:), Param(*), SolUf(:), Xel(:)
	integer, intent(in) :: NdimE, iFtype, constLaw, iFemType

	real*8 :: Fmu(NdimE,NdimE), dVng, PKaux(NdimE,NdimE) , GradUf(NdimE,NdimE)
	type(ptGaussClass) :: PtG
	integer :: NodG,NGP,pOrder,iBubble, iSimplex, nG
	
	PKmu = 0.0d0
	Vmu = 0.0d0
	
	call setFEMtype(iFemtype,NodG,pOrder,NGP,iSimplex,iBubble)
	
	call PtG%init(Xel,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	
	do nG = 1 , NGP
		call PtG%calcGradU(GradUf,SolUf,nG)
		
		call getF(Fmu,iFtype)
		
		Fmu = Fmu + GradUf
		
		PKaux = 0.0d0
		call calcPpk(PKaux,Fmu,NdimE,matPar,constLaw)			
!~ 		call damageModifyStress(PKaux,Param,nG)
		
		dVng=ptG%dV(nG)
		
		Vmu = Vmu + dVng
		PKmu = PKmu + dVng*PKaux
	end do	
	
	PKmu = PKmu/Vmu

end subroutine


subroutine calcGradU(GradU,Sol,dPhi_G,NdimE,NodElT,iDofSol)
	integer, intent(in) :: NdimE, NodElT, iDofSol 
	Real*8, intent(out) :: GradU(NdimE,NdimE)
	Real*8, intent(in) ::dPhi_G(NdimE,NodElT) , Sol(iDofSol*NodElt)
	integer :: i, j , e ,ep 
	
	GradU = 0.0d0
	do e = 1 , NodElT
		ep = (e-1)*iDofSol
		do i = 1, NdimE
			do j = 1, NdimE
				GradU(i,j) = GradU(i,j) + Sol(ep+i)*dPhi_G(j,e)
			end do
		end do
	end do	
end subroutine

subroutine calcU(U,Sol,Phi,NdimE,NodElT,iDofSol)
	integer, intent(in) :: NdimE,NodElT , iDofSol
	Real*8, intent(out) :: U(NdimE)
	Real*8, intent(in) :: Phi(NodElT) , Sol(iDofSol*NodElt)
	integer :: i, e ,ep 
	
	U= 0.0d0
	do e = 1 , NodElT
		ep = (e-1)*iDofSol
		do i = 1, NdimE
			U(i) = U(i) + Sol(ep+i)*Phi(e)
		end do
	end do	
	
end subroutine

real*8 function DGradGrad(DD,dPhiA,dPhiB,pA,pB,NdimE) result(rAux) !! correct order :: (D grad A) \cdot grad B
	real*8 , intent(in) :: DD(NdimE,NdimE,NdimE,NdimE),dPhiB(NdimE),dPhiA(NdimE)
	integer, intent(in) :: pA,pB, NdimE
	integer :: j,l
	
	rAux = 0.0d0
	
	do j=1,NdimE
	do l=1,NdimE
		rAux= rAux + DD(pB,j,pA,l)*dPhiB(j)*dPhiA(l) 
	end do
	end do
	
end function

real*8 function DGradPhiBdotGradPhiA(DD,dPhiB,dPhiA,p,q,NdimE) result(rAux) !! It would be more easy if it would q,p
	real*8 , intent(in) :: DD(NdimE,NdimE,NdimE,NdimE),dPhiB(NdimE),dPhiA(NdimE)
	integer, intent(in) :: p,q, NdimE
	integer :: j,l
	
	rAux = 0.0d0
	
	do j=1,NdimE
	do l=1,NdimE
		rAux= rAux + DD(p,j,q,l)*dPhiA(j)*dPhiB(l) 
	end do
	end do
	
end function
 
subroutine calcPpk(Ppk,F,NdimE,matPar,constLaw)
	integer , intent(in) :: NDimE,constLaw
	real*8, intent(in) :: MatPar(:), F(NDimE,NDimE)
	real*8, intent(out) :: Ppk(NDimE,NDimE)
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  Fp(NDimE,NDimE), inv2eps, energyPr , energyPl
	integer :: i,j 

	inv2eps = 0.5d0/eps
	
	Do i=1,NdimE  
	Do j=1,NdimE  
		
		Fp=F
		Fp(i,j) = Fp(i,j) + eps
		Call StrainEnergy(energyPr,Fp,MatPar,constLaw)

		Fp=F
		Fp(i,j) = Fp(i,j) - eps
		Call StrainEnergy(energyPl,Fp,MatPar,constLaw)

		Ppk(i,j) = (energyPr - energyPl)*inv2eps
	Enddo
	Enddo
	
end subroutine

subroutine calcD(D,F,NdimE,matPar,constLaw)
	integer , intent(in) :: NDimE,constLaw
	real*8, intent(in) :: MatPar(:), F(NDimE,NDimE)
	real*8, intent(out) :: D(NDimE,NDimE,NDimE,NDimE)
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  Fp(NDimE,NDimE), inveps2, inv4eps2, energy, energyPrr , energyPrl, energyPlr, energyPll
	integer :: i,j, k, l

	inv4eps2 = 0.25d0/(eps*eps)
	inveps2 = 1.0d0/(eps*eps)
	
	call strainEnergy(energy,F,MatPar,constLaw)
	
	Do i=1,NdimE  
	Do j=1,NdimE  
	Do k=1,NdimE  
	Do l=1,NdimE  
		
		if(i==k .and. j==l) then
			Fp=F
			Fp(i,j) = Fp(i,j) + eps
			Call StrainEnergy(energyPrr,Fp,MatPar,constLaw)
				
			Fp=F
			Fp(i,j) = Fp(i,j) - eps
			Call StrainEnergy(energyPll,Fp,MatPar,constLaw)
			
			D(i,j,i,j) = (energyPrr + energyPll - 2.0d0*energy)*inveps2
		else
			Fp=F
			Fp(i,j) = Fp(i,j) + eps
			Fp(k,l) = Fp(k,l) + eps
			Call StrainEnergy(energyPrr,Fp,MatPar,constLaw)
			
			Fp=F
			Fp(i,j) = Fp(i,j) + eps
			Fp(k,l) = Fp(k,l) - eps
			Call StrainEnergy(energyPrl,Fp,MatPar,constLaw)

			Fp=F
			Fp(i,j) = Fp(i,j) - eps
			Fp(k,l) = Fp(k,l) + eps
			Call StrainEnergy(energyPlr,Fp,MatPar,constLaw)
		
			Fp=F
			Fp(i,j) = Fp(i,j) - eps
			Fp(k,l) = Fp(k,l) - eps
			Call StrainEnergy(energyPll,Fp,MatPar,constLaw)
			
			D(i,j,k,l) = (energyPrr + energyPll - energyPrl - energyPlr)*inv4eps2
		end if
	end do
	end do
	end do
	end do
	
end subroutine

subroutine modifyToSpatial(Ds,Sigma,F,detF,NdimE)
	Real*8 , intent(out) :: sigma(NdimE,NdimE), Ds(NdimE,NdimE,NdimE,NdimE)
	Real*8 , intent(in) :: F(NdimE,NdimE),detF
	integer , intent(in) :: NDimE 
	
	Real*8 :: D(NdimE,NdimE,NdimE,NdimE), detFtemp
	integer :: i,j,k,l,m,n,p,q
	
	D = Ds
	
	sigma = (1.0d0/detF)*matmul(sigma,transpose(F))
	
	Ds = 0.0d0
	
	do i = 1,NdimE
	do j = 1,NdimE
	do k = 1,NdimE
	do l = 1,NdimE
	do p = 1,NdimE
	do q = 1,NdimE
		Ds(i,j,k,l) = Ds(i,j,k,l) +  D(i,p,k,q)*F(j,p)*F(l,q)/detF
	end do
	end do
	end do
	end do
	end do
	end do

end subroutine


subroutine calcFC(U,GradU,F,FT,FinvT,C,detF,Sol1,Phi,dPhi_G,NdimE,NodElT,iDofSol)
	integer, intent(in) :: NdimE,NodElT , iDofSol !!! Not exactly iDofT
	Real*8, intent(out) :: U(NdimE), GradU(NdimE,NdimE) , F(NdimE,NdimE), FT(NdimE,NdimE) , & 
								FinvT(NdimE,NdimE) , C(NdimE,NdimE) , detF
	Real*8, intent(in) ::Phi(NodElT), dPhi_G(NdimE,NodElT) , Sol1(iDofSol*NodElt)
	
	integer :: i, j , e ,ep 
	
	! Computes U and GradU
	U = 0.0d0
	GradU = 0.0d0
	do e = 1 , NodElT
		ep = (e-1)*iDofSol
		do i = 1, NdimE
			U(i) = U(i) + Sol1(ep+i)*Phi(e)
			do j = 1, NdimE
!~ 				if(i==3) write(*,*) Sol1(ep+i), dPhi_G(j,e)
				GradU(i,j) = GradU(i,j) + Sol1(ep+i)*dPhi_G(j,e)
			end do
		end do
	end do	
	
!~ 	call printVec(Sol1(1:9))
!~ 	call printMat(dPhi_G)
!~ 	! Computes F = I + GradU + F0 
	call addEye(F) !! F = F0 + I
	F = F + GradU !!! F might be initialized before as zeros or with a given deformation
	FT = transpose(F) !~ 	! Computes  FT = (F)^T
	C = matmul(FT,F) ! Cauchy-Grenn tensor
	call MatInv(FinvT,detF,FT)  ! FinvT= (FT)^-1 , defF = det(FT) = det(F)
	
	return
end subroutine


subroutine VirtualPower()
	use ptsGaussLib2
	
!~ 	integer , parameter ::  NdimE = 2 , iFemType = 2 !! For linear quadrangles
	integer , parameter ::  NdimE = 3 , iFemType = 6 !! tetrahedra
	integer :: pOrder , NGP, NodG, iSimplex, iBubble, nG, constLaw
	real*8 :: MatPar(3), alpha
	real*8 , allocatable :: Xel(:), SolU(:), Xs(:), Xm(:), SolV(:)
	
	Real*8 :: dVs , detF0s, detFs 
	Real*8 , dimension(NdimE,NdimE) :: GradUs , GradU0s, GradVs, DsGradUs, sigma
	Real*8 , dimension(NdimE,NdimE) :: Fs, F0s, Fbars,  tenQsGradUs, tenQsGradU0s
	Real*8 , dimension(NdimE,NdimE,NdimE,NdimE) :: Ds , tenQs
	type(ptGaussClass) :: PtGs, PtG0s
	real*8 :: PotAs, PotQs, PotQ0s, PotTs 

	Real*8 :: dVm , detF0m, detFm 
	Real*8 , dimension(NdimE,NdimE) :: GradUm , GradU0m, GradVm, DmGradUm, Ppk
	Real*8 , dimension(NdimE,NdimE) :: Fm, F0m, Fbarm,  tenQmGradUm, tenQmGradU0m
	Real*8 , dimension(NdimE,NdimE,NdimE,NdimE) :: Dm , tenQm, tenQ0m
	type(ptGaussClass) :: PtGm, PtG0m
	real*8 :: PotAm, PotQm, PotQ0m, PotTm 

	constLaw = 4
	MatPar = (/100.0d0, 200.0d0, 5000.0d0/)
	
	if(NdimE == 2) then
		alpha = 0.5d0
	else if(NdimE == 3) then
		alpha = 1.0d0/3.0d0
	end if
	
	call setFEMtype(iFEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
	write(0,*) NGP
	
	allocate(Xel(NdimE*NodG),SolU(NdimE*NodG),Xs(NdimE*NodG),Xm(NdimE*NodG),SolV(NdimE*NodG))
	
	!! quadrangles
!~ 	Xel = (/-1.0d0,-1.0d0,1.0d0,-1.0d0,1.0d0, 1.0d0,-1.0d0, 1.0d0/)
!~ 	SolU = (/0.01d0,0.02d0,0.03d0,0.04d0,0.05d0, 0.06d0,0.07d0, 0.08d0/)
!~ 	SolV = (/-1.0d0,1.0d0,2.0d0,3.0d0,-0.1d0, 4.0d0,5.0d0, -5.0d0/)
!~ 	
	!! tetrahedra !!! linear, is not well posed the method
	Xel  = (/-0.02d0, -0.02d0, -0.02d0, 1.00d0, 0.00d0, 0.00d0, 0.00d0, 1.00d0, 0.00d0, 0.00d0, 0.00d0, 1.00d0 /)
	SolU = (/0.02d0, 0.02d0, 0.03d0, 0.04d0, 0.05d0, 0.06d0, 0.07d0, 0.08d0, 0.02d0, 0.03d0, 0.04d0, -0.00d0 /)
	SolV = (/0.01d0, 0.03d0, -0.03d0, -0.04d0, 0.08d0, 0.01d0, -0.07d0, -0.09d0, -0.01d0, 0.02d0, 0.03d0, -0.05d0 /)
	
	Xs = Xel + SolU
	Xm = Xel
	
	call PtGs%init(Xs,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	call PtG0s%init(Xs,NodG,NdimE,1,pOrder,iBubble,iSimplex)
	call PtGm%init(Xm,NodG,NdimE,NGP,pOrder,iBubble, iSimplex)
	call PtG0m%init(Xm,NodG,NdimE,1,pOrder,iBubble,iSimplex)
	
	call numprint(ptGs%dV)
	call numprint(ptGm%dV)
	
	
	nG = 1
	call PtG0s%calcGradU(GradU0s,SolU,nG)
	F0s = 0.0d0
	call computeFs(F0s,detF0s,gradU0s,NdimE)

	nG = 1
	call PtG0m%calcGradU(GradU0m,SolU,nG)
	F0m = deltaKron(1:NdimE,1:NdimE) + GradU0m
	call calcI3(detF0m,F0m)
	
	PotAs = 0.0d0; PotQs = 0.0d0; PotQ0s = 0.0d0; PotTs = 0.0d0
	PotAm = 0.0d0; PotQm = 0.0d0; PotQ0m = 0.0d0; PotTm = 0.0d0 
			
	Do nG = 1, NGP ! LoopGauss
	
		call PtGs%calcGradU(GradUs,SolU,nG)
		call PtGs%calcGradU(GradVs,SolV,nG)
		call PtGm%calcGradU(GradUm,SolU,nG)
		call PtGm%calcGradU(GradVm,SolV,nG)
		
		Fs = 0.0d0
		call computeFs(Fs,detFs,gradUs,NdimE)
		
		Fm = deltaKron(1:NdimE,1:NdimE) + GradUm
		call calcI3(detFm,Fm)
		
		write(0,*) detF0s/detFs
		
		Fbars = ((detF0s/detFs)**alpha)*Fs
		Fbarm = ((detF0m/detFm)**alpha)*Fm
		
		call calcPpk(sigma,Fbars,NdimE,matPar,constLaw)		
		call calcD(Ds,Fbars,NdimE,matPar,constLaw)
		call modifyToSpatial(Ds,sigma,Fbars,detF0s,NdimE)
		
		call calcPpk(Ppk,Fbarm,NdimE,matPar,constLaw)		
		call calcD(Dm,Fbarm,NdimE,matPar,constLaw)
			
		call build_tenQs(Ds,sigma,tenQs,NdimE)
		call buildFbarTensors(Dm,Ppk,tenQm,tenQ0m,Fm,F0m,NdimE)

		dVs=ptGs%dV(nG)
		dVm=ptGm%dV(nG)
		
		call T4xT2(DsGradUs,Ds,GradUs)
		call T4xT2(tenQsGradUs,tenQs,GradUs)
		call T4xT2(tenQsGradU0s,tenQs,GradU0s)

		call T4xT2(DmGradUm,Dm,GradUm)
		call T4xT2(tenQmGradUm,tenQm,GradUm)
		call T4xT2(tenQmGradU0m,tenQ0m,GradU0m)
		
		PotAs = PotAs + dot_product2(DsGradUs,GradVs)*dvs
		PotQs = PotQs + dot_product2(tenQsGradUs,GradVs)*dvs
		PotQ0s = PotQ0s + dot_product2(tenQsGradU0s,GradVs)*dvs
		PotTs = PotTs + dot_product2(sigma,GradVs)*dvs
	
		PotAm = PotAm + dot_product2(DmGradUm,GradVm)*dvm
		PotQm = PotQm + dot_product2(tenQmGradUm,GradVm)*dvm
		PotQ0m = PotQ0m + dot_product2(tenQmGradU0m,GradVm)*dvm	
		PotTm = PotTm + dot_product2(Ppk,GradVm)*dvm
	
		
	end do
		
	write(0,*) PotAs, PotQs, PotQ0s, PotTs
	write(0,*) PotAm, PotQm, PotQ0m, PotTm
	pause
end subroutine

subroutine buildFbarTensors(Astar,PpkStar,Qstar,Q0star,F,F0,NdimE)
	
	real*8 , intent(inout) :: Astar(NdimE,NdimE,NdimE,NdimE), PpkStar(NdimE,NdimE)
	real*8 , intent(inout) :: Qstar(NdimE,NdimE,NdimE,NdimE), Q0star(NdimE,NdimE,NdimE,NdimE)
	real*8 , intent(in) :: F(NdimE,NdimE), F0(NdimE,NdimE)
	integer , intent(in) :: NdimE
	
	real*8 :: Finv(NdimE,NdimE), F0inv(NdimE,NdimE), Jratio, detF, detF0, alpha , beta, gamma, lamb
	integer :: i,j,k,l,p,q	

	call MatInv(Finv,detF,F) !! maybe detF0 is not correct because is 2D
	call MatInv(F0inv,detF0,F0) !! maybe detF0 is not correct because is 2D
	
	if(NdimE == 2) then
		alpha = -0.5d0
		beta = 0.0d0
		gamma = 0.5d0
		lamb = -0.5d0 
	else if(NdimE == 3) then
		alpha =  -2.0d0/3.0d0
		beta = -1.0d0/3.0d0
		gamma = 1.0d0/3.0d0
		lamb = -2.0d0/3.0d0
	end if
	

	Jratio = detF0/detF
	
	PpkStar = (Jratio**alpha)*PpkStar
	Astar = (Jratio**beta)*Astar
	
	Qstar = 0.0d0
	Q0star = 0.0d0
	
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		do p = 1 , NdimE
		do q = 1 , NdimE
			Qstar(i,j,k,l) = Qstar(i,j,k,l) + gamma*Astar(i,j,p,q)*F(p,q)*Finv(l,k) !! for 2D
			Q0star(i,j,k,l) = Q0star(i,j,k,l) + gamma*Astar(i,j,p,q)*F(p,q)*F0inv(l,k) !! for 2D
		end do
		end do
		Qstar(i,j,k,l) = Qstar(i,j,k,l) + lamb*PpkStar(i,j)*Finv(l,k) 
		Q0star(i,j,k,l) = Q0star(i,j,k,l) + lamb*PpkStar(i,j)*F0inv(l,k) 
	end do
	end do
	end do
	end do
	
end subroutine

subroutine build_tenQs(Ds,sigma,tenQs,NdimE)	
	real*8 , intent(inout) :: Ds(NdimE,NdimE,NdimE,NdimE), sigma(NdimE,NdimE)
	real*8 , intent(inout) :: tenQs(NdimE,NdimE,NdimE,NdimE)
	integer , intent(in) :: NdimE
	
	integer :: i,j,k,l,p,q
	real*8 :: gamma, lamb
	
	tenQs = 0.0d0
	
	if(NdimE == 2) then
		gamma = 0.5d0
		lamb = -0.5d0 
	else if(NdimE == 3) then
		gamma = 1.0d0/3.0d0
		lamb = -2.0d0/3.0d0
	end if
	
	
	do i = 1 , NdimE
	do j = 1 , NdimE
	do k = 1 , NdimE
	do l = 1 , NdimE
		do p = 1 , NdimE
		do q = 1 , NdimE
			tenQs(i,j,k,l) = tenQs(i,j,k,l) + gamma*Ds(i,j,p,q)*deltaKron(p,q)*deltaKron(k,l)
		end do
		end do
		tenQs(i,j,k,l) = tenQs(i,j,k,l) + lamb*sigma(i,j)*deltaKron(k,l) 
	end do
	end do
	end do
	end do
	
end subroutine


end module

