module linElasLib
use funcAux
use damageNewLib
implicit none

public buildElasticityTensor, buildElasticityTensorMatSym, getEpsMatSym, testElasticity, strainEnergyLinElas

contains

subroutine buildElasticityTensor(C,E,nu,NdimE)
	Real*8 , intent(out) :: C(NdimE,NdimE,NdimE,NdimE) 
	Real*8 , intent(in) :: E, nu
	integer , intent(in) :: NdimE
	integer :: i,j,k,l
	
	Real*8 :: lamb, mu
	
	lamb = E*nu/((1.0d0+nu)*(1.0d0 - 2.0d0*nu))
	mu = E/(2.0d0*(1.0d0 + nu))
	
	do i = 1,NdimE
	do j = 1,NdimE
	do k = 1,NdimE
	do l = 1,NdimE
	
		C(i,j,k,l) = 2.0d0*mu*deltaKron(i,k)*deltaKron(j,l) + lamb*deltaKron(i,j)*deltaKron(k,l)
	
	end do
	end do
	end do
	end do

end subroutine

subroutine buildElasticityTensorMatSym(D,E,nu,NdimE) ! Following Hughes pág 83
	Real*8 , intent(out) :: D(NdimE*NdimE-1,NdimE*NdimE-1) 
	Real*8 , intent(in) :: E, nu
	integer , intent(in) :: NdimE
	integer :: i,j,k,l
	
	Real*8 :: lamb, mu
	
	lamb = E*nu/((1.0d0+nu)*(1.0d0 - 2.0d0*nu))
	mu = E/(2.0d0*(1.0d0 + nu))
	
!~ 	lamb = 2.0d0*lamb*mu/(lamb + 2.0d0*mu)
	
	D = 0.0d0
	
	D(1,1) = lamb + 2.0d0*mu
	D(2,2) = lamb + 2.0d0*mu
	D(3,3) = mu
	D(1,2) = lamb
	D(2,1) = lamb

end subroutine

subroutine getEpsMatSym(epsMat,GradU)
	real*8, intent(out) :: epsMat(:)
	real*8, intent(in) :: GradU(:,:)
	
	epsMat(1) = GradU(1,1)
	epsMat(2) = GradU(2,2)
	epsMat(3) = GradU(1,2) + GradU(2,1)
end subroutine

subroutine testElasticity()
	use funcAux 
	implicit none
	
	integer , parameter :: NdimE = 2, NdimEMat = 3 
	real*8 :: E,nu, D(NdimE,NdimE,NdimE,NdimE),  sigma(NdimE,NdimE), eps(NdimE,NdimE), &
			  Dmat(NdimEMat,NdimEMat), epsMat(NdimEMat) ,  sigmaMat(NdimEMat), &
			  Dmat2(NdimE*NdimE,NdimE*NdimE), epsMat2(NdimE*NdimE) ,  sigmaMat2(NdimE*NdimE)
	
	
	eps(1,1) = 0.1d0
	eps(1,2) = -0.1d0
	eps(2,1) = eps(1,2)
	eps(2,2) = 0.1d0
	
	epsMat(1) = eps(1,1)
	epsMat(2) = eps(2,2)
	epsMat(3) = 2.0d0*eps(1,2)
	
	E = 100.0d0
	nu = 0.3d0
	
	call voigtTen2toVec(epsMat2,eps)
	call buildElasticityTensor(D,E,nu,NdimE)
 	call voigtTen4toTen2(Dmat2,D)
 	call buildElasticityTensorMatSym(Dmat,E,nu,NdimE)
	
	call T4xT2(sigma,D,eps)
	sigmaMat = matmul(Dmat,epsMat)
	sigmaMat2 = matmul(Dmat2,epsMat2)
	
	call numprint(sigma)
	call numprint(sigmaMat)
	call numprint(sigmaMat2)

end subroutine

subroutine strainEnergyLinElas(energy,GradU,matpar)
	real*8, intent(in) :: GradU(:,:), matpar(:)
	real*8 , intent(out) :: energy
	
	integer, parameter :: NdimE = 2
	real*8 :: E, nu, Dmat(NdimE*NdimE-1,NdimE*NdimE-1), epsMat(NdimE*NdimE - 1)
	
	E = matpar(1)
	nu = matpar(2)
	
	call getEpsMatSym(epsMat,GradU)
	call buildElasticityTensorMatSym(Dmat,E,nu,NdimE)
	
	energy = 0.5d0*dot_product(epsMat,matmul(Dmat,epsMat))
	
end subroutine 

subroutine damageModifyTangentLinElas(Dmat,Dmat0,damagePar,Param, nG)
	real*8, intent(inout) :: Dmat(:,:)
	real*8 , intent(in) :: Dmat0(:,:), damagePar(:), Param(*)
	integer , intent(in) :: nG
	type(damageState) :: sd

	call sd%loadDamage(Param, nG, .False.)
	
!~ 	write(0,*) "DAMAGE ==============>" , sd%damage

	Dmat = (1.0d0-sd%damage)*Dmat0 

end subroutine

end module
