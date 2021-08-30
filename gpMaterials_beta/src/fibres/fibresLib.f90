module fibresLib

use funcAux
use fibresModelsLib

implicit none

public calcSfib, calcDSfib, setFibreKinematic, setFibreKinematicSimple, getPKmuFromfib

!~ 	call calcSfib(Sfib,qfib,matPar,constLaw)		
!~ 	call calcDSfib(DSfib,qfib,matPar,constLaw)	

contains

subroutine getPKmuFromfib(PKfib,Vfib,SolUf,Xel,matpar,Param,constLaw,iFtype,NdimE)

	use damageFibLib
	
	real*8 , intent(out) :: PKfib(:,:) , Vfib
	real*8, intent(in) :: matPar(:), Param(*), SolUf(:), Xel(:)
	integer, intent(in) :: NdimE, iFtype, constLaw
	Real*8 :: Sfib(NdimE), AreaFib, DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib
	
	AreaFib = matpar(3)
	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolUf,Xel,iFtype)
	call calcSfib(Sfib,qfib,matPar,constLaw,NDimE)
!~ 	call damageModifySfib(Sfib, Param)
	Vfib = AreaFib*Lfib
		
	PKfib(1,1) = Vfib*Sfib(1)*afib(1)
	PKfib(1,2) = Vfib*Sfib(1)*afib(2)
	PKfib(2,1) = Vfib*Sfib(2)*afib(1)
	PKfib(2,2) = Vfib*Sfib(2)*afib(2)	

end subroutine

real*8 function getLambdaFibre(SolU,Xel,iFtype,NdimE) result(lambda) 
	Real*8 , intent(in) :: SolU(:) , Xel(:)
	integer , intent(in) :: NdimE, iFtype
	Real*8  :: DeltaU(NdimE), afib(NdimE), qfib(NdimE), Lfib
	
	call setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)
	
	lambda = norm2(qfib)
	
end function

subroutine get_afib(afib,Lfib,Xel)
	Real*8 , intent(out) :: afib(:), Lfib
	Real*8 , intent(in) :: Xel(:)
			
	afib(1) = Xel(3) - Xel(1)
	afib(2) = Xel(4) - Xel(2)  
	
	Lfib = norm2(afib)
	afib = afib/Lfib

end subroutine

subroutine get_afibbar(afibbar,Xel)
	Real*8 , intent(out) :: afibbar(:)
	Real*8 , intent(in) :: Xel(:)
			
	afibbar(1) = 0.5d0*(Xel(3) + Xel(1))
	afibbar(2) = 0.5d0*(Xel(4) + Xel(2))  
	
end subroutine


subroutine setFibreKinematic(qfib,DeltaU,afib,Lfib,SolU,Xel,iFtype)
	use globalVariables, only : getF, NdimE
	
	Real*8 , intent(out) :: DeltaU(:), afib(:), qfib(:), Lfib
	Real*8 , intent(in) :: SolU(:) , Xel(:)
	integer , intent(in) :: iFtype
	real*8 :: F(NdimE,NdimE)
	
	call getF(F,iFtype)

	call get_afib(afib,Lfib,Xel)
	
	DeltaU(1) = SolU(3) - SolU(1)
	DeltaU(2) = SolU(4) - SolU(2)
	
	qfib = matmul(F,afib) + DeltaU/Lfib	

end subroutine

subroutine setFibreKinematicSimple(qfib,DeltaU,afib,Lfib,SolU,iFtype)
	use globalVariables, only : getF, NdimE
	
	Real*8 , intent(out) :: DeltaU(:), qfib(:)
	Real*8 , intent(in) :: SolU(:) , afib(:), Lfib
	integer , intent(in) :: iFtype
	real*8 :: F(NdimE,NdimE)
	
	call getF(F,iFtype)
!~ 	call numprint(F)
	
	
	DeltaU(1) = SolU(3) - SolU(1)
	DeltaU(2) = SolU(4) - SolU(2)
	
	qfib = matmul(F,afib) + DeltaU/Lfib	

end subroutine

subroutine calcSfib(Sfib,qfib,matPar,constLaw,NdimE)
	integer , intent(in) :: constLaw, NDimE
	real*8, intent(in) :: MatPar(:), qfib(NDimE)
	real*8, intent(out) :: Sfib(NDimE)
	logical , parameter :: isAnalitic = .true.
	
	if(isAnalitic) then
		call SfibAnalytic(Sfib,qfib,matPar,constLaw)
	else
		call calcSfibNumeric(Sfib,qfib,matPar,constLaw,NdimE)
	end if

end subroutine


subroutine calcDSfib(DSfib,qfib,matPar,constLaw,NdimE)
	integer , intent(in) :: constLaw, NDimE
	real*8, intent(in) :: MatPar(:), qfib(NDimE)
	real*8, intent(out) :: DSfib(NDimE,NdimE)
	logical , parameter :: isAnalitic = .false.
	
	if(isAnalitic) then
		call DSfibAnalytic(DSfib,qfib,matPar,constLaw)
	else
		call calcDSfibNumeric(DSfib,qfib,matPar,constLaw,NdimE)
	end if

end subroutine

subroutine calcSx(Sx,Ex,matPar,constLaw)
	integer , intent(in) :: constLaw
	real*8, intent(in) :: MatPar(:), Ex
	real*8, intent(out) :: Sx
	logical , parameter :: isAnalitic = .false.
	
	if(isAnalitic) then
		call SxAnalytic(Sx,Ex,matPar,constLaw)
	else
		call calcSxNumeric(Sx,Ex,matPar,constLaw)
	end if

end subroutine


subroutine calcDSx(DSx,Ex,matPar,constLaw)
	integer , intent(in) :: constLaw
	real*8, intent(in) :: MatPar(:), Ex
	real*8, intent(out) :: DSx
	logical , parameter :: isAnalitic = .false.
	
	if(isAnalitic) then
		call DSxAnalytic(DSx,Ex,matPar,constLaw)
	else
		call calcDSxNumeric(DSx,Ex,matPar,constLaw)
	end if

end subroutine


subroutine calcSfibNumeric(Sfib,qfib,matPar,constLaw,NdimE)
	integer , intent(in) :: constLaw, NDimE
	real*8, intent(in) :: MatPar(:), qfib(NDimE)
	real*8, intent(out) :: Sfib(NDimE)
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  qfibP(NDimE), inv2eps, energyPr , energyPl
	integer :: i 

	inv2eps = 0.5d0/eps
	
	Do i=1,NdimE    
		qfibP=qfib
		qfibP(i) = qfibP(i) + eps
		Call strainEnergy(energyPr,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) - eps
		Call strainEnergy(energyPl,qfibP,MatPar,constLaw)

		Sfib(i) = (energyPr - energyPl)*inv2eps
	Enddo
	
end subroutine

subroutine calcDSfibNumeric(DSfib,qfib,matPar,constLaw,NdimE)
	integer , intent(in) :: constLaw, NDimE
	real*8, intent(in) :: MatPar(:), qfib(NDimE)
	real*8, intent(out) :: DSfib(NDimE,NdimE)
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  qfibP(NDimE), inveps2, inv4eps2, energy, energyPr , energyPl,  &
				energyPrr , energyPrl , energyPlr, energyPll 
	integer :: i,j 

	Call strainEnergy(energy,qfib,MatPar,constLaw)

	inveps2 = 1.0d0/(eps*eps)
	inv4eps2 = 0.25d0/(eps*eps)	

	Do i=1,NdimE    
		qfibP=qfib
		qfibP(i) = qfibP(i) + eps
		Call strainEnergy(energyPr,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) - eps
		Call strainEnergy(energyPl,qfibP,MatPar,constLaw)

		DSfib(i,i) = (energyPr + energyPl - 2.0d0*energy)*inveps2
	Enddo

	Do i=1,NdimE    
	Do j=i+1,NdimE        
		qfibP=qfib
		qfibP(i) = qfibP(i) + eps
		qfibP(j) = qfibP(j) + eps
		Call strainEnergy(energyPrr,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) + eps
		qfibP(j) = qfibP(j) - eps
		Call strainEnergy(energyPrl,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) - eps
		qfibP(j) = qfibP(j) + eps
		Call strainEnergy(energyPlr,qfibP,MatPar,constLaw)

		qfibP=qfib
		qfibP(i) = qfibP(i) - eps
		qfibP(j) = qfibP(j) - eps
		Call strainEnergy(energyPll,qfibP,MatPar,constLaw)

		DSfib(i,j) = (energyPrr + energyPll - energyPrl - energyPlr)*inv4eps2
		DSfib(j,i) = DSfib(i,j)
	Enddo
	Enddo
 	
end subroutine

subroutine calcSxNumeric(Sx,Ex,matPar,constLaw)
	integer , intent(in) :: constLaw
	real*8, intent(in) :: MatPar(:), Ex
	real*8, intent(out) :: Sx
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  Ex_p, inv2eps, energyPr , energyPl

	inv2eps = 0.5d0/eps
	   
	Ex_p=Ex + eps
	Call strainEnergyEx(energyPr,Ex_p,MatPar,constLaw)

	Ex_p=Ex - eps	
	Call strainEnergyEx(energyPl,Ex_p,MatPar,constLaw)

	Sx = (energyPr - energyPl)*inv2eps

end subroutine

subroutine calcDSxNumeric(DSx,Ex,matPar,constLaw)
	integer , intent(in) :: constLaw
	real*8, intent(in) :: MatPar(:), Ex
	real*8, intent(out) :: DSx
	real*8 , parameter :: eps = 1.d-4
	real*8 ::  Ex_p, inveps2, energy, energyPr , energyPl
	integer :: i,j 

	Call strainEnergyEx(energy,Ex,MatPar,constLaw)

	inveps2 = 1.0d0/(eps*eps)
	
	Ex_p=Ex + eps
	Call strainEnergyEx(energyPr,Ex_p,MatPar,constLaw)

	Ex_p=Ex - eps	
	Call strainEnergyEx(energyPl,Ex_p,MatPar,constLaw)

	DSx = (energyPr + energyPl - 2.0d0*energy)*inveps2

end subroutine

!~ subroutine calcGradU(GradU,Sol,dPhi_G,NdimE,NodElT,iDofSol)
!~ 	integer, intent(in) :: NdimE, NodElT, iDofSol 
!~ 	Real*8, intent(out) :: GradU(NdimE,NdimE)
!~ 	Real*8, intent(in) ::dPhi_G(NdimE,NodElT) , Sol(iDofSol*NodElt)
!~ 	integer :: i, j , e ,ep 
!~ 	
!~ 	GradU = 0.0d0
!~ 	do e = 1 , NodElT
!~ 		ep = (e-1)*iDofSol
!~ 		do i = 1, NdimE
!~ 			do j = 1, NdimE
!~ 				GradU(i,j) = GradU(i,j) + Sol(ep+i)*dPhi_G(j,e)
!~ 			end do
!~ 		end do
!~ 	end do	
!~ end subroutine
!~ 
!~ subroutine calcU(U,Sol,Phi,NdimE,NodElT,iDofSol)
!~ 	integer, intent(in) :: NdimE,NodElT , iDofSol
!~ 	Real*8, intent(out) :: U(NdimE)
!~ 	Real*8, intent(in) :: Phi(NodElT) , Sol(iDofSol*NodElt)
!~ 	integer :: i, e ,ep 
!~ 	
!~ 	U= 0.0d0
!~ 	do e = 1 , NodElT
!~ 		ep = (e-1)*iDofSol
!~ 		do i = 1, NdimE
!~ 			U(i) = U(i) + Sol(ep+i)*Phi(e)
!~ 		end do
!~ 	end do	
!~ 	
!~ end subroutine
 


!~ subroutine calcD(D,F,NdimE,matPar,constLaw)
!~ 	integer , intent(in) :: NDimE,constLaw
!~ 	real*8, intent(in) :: MatPar(:), F(NDimE,NDimE)
!~ 	real*8, intent(out) :: D(NDimE,NDimE,NDimE,NDimE)
!~ 	real*8 , parameter :: eps = 1.d-4
!~ 	real*8 ::  Fp(NDimE,NDimE), inveps2, inv4eps2, energy, energyPrr , energyPrl, energyPlr, energyPll
!~ 	integer :: i,j, k, l
!~ 
!~ 	inv4eps2 = 0.25d0/(eps*eps)
!~ 	inveps2 = 1.0d0/(eps*eps)
!~ 	
!~ 	call strainEnergy(energy,F,MatPar,constLaw)
!~ 	
!~ 	Do i=1,NdimE  
!~ 	Do j=1,NdimE  
!~ 	Do k=1,NdimE  
!~ 	Do l=1,NdimE  
!~ 		
!~ 		if(i==k .and. j==l) then
!~ 			Fp=F
!~ 			Fp(i,j) = Fp(i,j) + eps
!~ 			Call StrainEnergy(energyPrr,Fp,MatPar,constLaw)
!~ 				
!~ 			Fp=F
!~ 			Fp(i,j) = Fp(i,j) - eps
!~ 			Call StrainEnergy(energyPll,Fp,MatPar,constLaw)
!~ 			
!~ 			D(i,j,i,j) = (energyPrr + energyPll - 2.0d0*energy)*inveps2
!~ 		else
!~ 			Fp=F
!~ 			Fp(i,j) = Fp(i,j) + eps
!~ 			Fp(k,l) = Fp(k,l) + eps
!~ 			Call StrainEnergy(energyPrr,Fp,MatPar,constLaw)
!~ 			
!~ 			Fp=F
!~ 			Fp(i,j) = Fp(i,j) + eps
!~ 			Fp(k,l) = Fp(k,l) - eps
!~ 			Call StrainEnergy(energyPrl,Fp,MatPar,constLaw)
!~ 
!~ 			Fp=F
!~ 			Fp(i,j) = Fp(i,j) - eps
!~ 			Fp(k,l) = Fp(k,l) + eps
!~ 			Call StrainEnergy(energyPlr,Fp,MatPar,constLaw)
!~ 		
!~ 			Fp=F
!~ 			Fp(i,j) = Fp(i,j) - eps
!~ 			Fp(k,l) = Fp(k,l) - eps
!~ 			Call StrainEnergy(energyPll,Fp,MatPar,constLaw)
!~ 			
!~ 			D(i,j,k,l) = (energyPrr + energyPll - energyPrl - energyPlr)*inv4eps2
!~ 		end if
!~ 	end do
!~ 	end do
!~ 	end do
!~ 	end do
!~ 	
!~ end subroutine



end module

