module damageCriteria
use funcAux
implicit none

real*8 :: Ikmax , J_k , J_km , dJrel , Idam, measOmega

contains

subroutine sameAsDissertation2(damageNew,AlphaNew,XLL,Sol1,EvolPar,matPar,damageOld,alphaOld,nDim,&
							constLaw,NodElT,iDofT,iShiftV,iShiftW,Param,lengthParam)
							
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftW,lengthParam  ! all integers
    Real*8, intent(in) :: XLL(*),Sol1(*),EvolPar(:),matPar(:), damageOld, alphaOld,Param(*)
    Real*8, intent(out) :: damageNew,AlphaNew

    Real*8 :: F, vardamage, deltaDamage, dampFac, ScaleFac , W(Ndim) , Ik , &
				Ikrel, threshold, gamma, kappa, dhat, avgDam, beta
    Real*8 , parameter :: MaxDamage = 0.965d0
    Real*8 :: maxDeltaDamage , minDeltaDamage , maxDamage_k
    
    
    maxDeltaDamage = 0.025d0
    minDeltaDamage = -0.025d0
    
    gamma = EvolPar(1)
    kappa = EvolPar(2)
    threshold = EvolPar(3)
    
	Ik = Param(lengthParam + 6)
	Ikrel = Ik/Ikmax
	
	avgDam = Idam/measOmega
	dJrel = (J_km - J_k)/J_k
		
	DeltaDamage = kappa*J_k*Ik

	write(0,*) "DeltaDamageNew = " , deltaDamage ,  "J_k= " , J_k , "Idam= " , avgDam, "dJrel= " , dJRel
	
	if(deltaDamage > maxDeltaDamage) then 
		deltaDamage = maxDeltaDamage
	else if(deltaDamage < minDeltaDamage) then
		deltaDamage = minDeltaDamage
	end if
	
	DamageNew = DamageOld + deltaDamage
		
	if(DamageNew>MaxDamage)  DamageNew = maxDamage
	if(DamageNew<0.0)  DamageNew = 0.0	
	
	AlphaNew = 0.0d0 !! just for compatibility with ifort
	       
end subroutine

subroutine sameAsDissertation(damageNew,AlphaNew,XLL,Sol1,EvolPar,matPar,damageOld,alphaOld,nDim,&
							constLaw,NodElT,iDofT,iShiftV,iShiftW,Param,lengthParam)
							
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftW,lengthParam  ! all integers
    Real*8, intent(in) :: XLL(*),Sol1(*),EvolPar(:),matPar(:), damageOld, alphaOld,Param(*)
    Real*8, intent(out) :: damageNew,AlphaNew

    Real*8 :: F, vardamage, deltaDamage, dampFac, ScaleFac , W(Ndim) , Ik , &
				Ikrel, threshold, gamma, kappa, dhat, avgDam, beta
    Real*8 , parameter :: MaxDamage = 0.965d0
    Real*8 :: maxDeltaDamage , minDeltaDamage , maxDamage_k
    
    
    maxDeltaDamage = 0.1d0
    minDeltaDamage = -0.1d0
    
    gamma = EvolPar(1)
    kappa = EvolPar(2)
    threshold = EvolPar(3)
    
	Ik = Param(lengthParam + 6)
	Ikrel = Ik/Ikmax
	
	avgDam = Idam/measOmega
	dJrel = (J_km - J_k)/J_k
	
!~ 	beta = -10.0d0 * dlog(0.1d0/maxDamage)
!~ 	maxDeltaDamage = maxDamage*dexp(-beta*avgDam)

	if(Ik>threshold) then
		dhat = gamma*(1.0d0-damageOld)
	else if(Ik<0.0d0) then
		dhat = -gamma*damageOld
	else
		dhat = 0.0d0
	end if
		
!~ 	maxDamage_k = 0.1d0 + 5.0d0*avgDam
!~ 	if(maxDamage_k>MaxDamage) maxDamage_k = MaxDamage 
	
!~ 	kappa = 1.5d0*maxDamage_k/Ikmax
	DeltaDamage = kappa*J_k*dhat

	write(0,*) "DeltaDamageNew = " , deltaDamage ,  "J_k= " , J_k , "Idam= " , avgDam, "dJrel= " , dJRel
	
	if(deltaDamage > maxDeltaDamage) then 
		deltaDamage = maxDeltaDamage
	else if(deltaDamage < minDeltaDamage) then
		deltaDamage = minDeltaDamage
	end if
	
	DamageNew = DamageOld + deltaDamage
		
	if(DamageNew>MaxDamage)  DamageNew = maxDamage
	if(DamageNew<0.0)  DamageNew = 0.0	
	
	AlphaNew = 0.0d0 !! just for compatibility with ifort
	       
end subroutine


subroutine Feijoo(damageNew,AlphaNew,XLL,Sol1,EvolPar,matPar,damageOld,alphaOld,nDim,&
							constLaw,NodElT,iDofT,iShiftV,iShiftW,Param,lengthParam)
							
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftW,lengthParam  ! all integers
    Real*8, intent(in) :: XLL(*),Sol1(*),EvolPar(:),matPar(:), damageOld, alphaOld,Param(*)
    Real*8, intent(out) :: damageNew,AlphaNew

    Real*8 :: F, vardamage, deltaDamage, dampFac, ScaleFac , W(Ndim) , sens, sensRel , threshold
    Real*8 , parameter :: MaxDamage = 0.98d0 ,  maxDeltaDamage = 0.25d0
    
    dampFac = EvolPar(1)
    ScaleFac = EvolPar(2)
    threshold = EvolPar(3)
    
	sens = Param(lengthParam + 6)
!~ 	sensRel = dsign(dsqrt(dabs(sens/MaxSens)),sens)
	
	
	if(sens>0.0d0) then
		sensRel = dsqrt(sens*(1.0d0-damageOld)/Ikmax)
	else
		sensRel = -dsqrt(-damageOld*sens/Ikmax)
	end if
	
!~ 	sensRel = sens/MaxSens
	deltaDamage = dsqrt(J_k)*sensRel*maxDeltaDamage
	
		
!~ 	deltaDamage = (1.0d0 - 3.0d0*Jcost)*ScaleFac*sens
!~ 	deltaDamage = Jcost*ScaleFac*sensRel

!~ 	if(sensRel>threshold .or. sensRel<-threshold) then
!~ 		deltaDamage = Jcost*ScaleFac*sensRel
!~ 	else
!~ 		deltaDamage = 0.0d0
!~ 	end if
	
!~     if(sensRel>threshold) then
!~ 		deltaDamage = ScaleFac*(1.0d0 - damageOld)
!~ 	else if(sensRel<-threshold) then
!~ 		deltaDamage = -ScaleFac*damageOld
!~ 	else
!~ 		deltaDamage = 0.0d0
!~ 	end if
!~     
    write(0,*) "DeltaDamageNew = " , deltaDamage ,  "costFunc= " , J_k
	
	if(deltaDamage > maxDeltaDamage) then 
		deltaDamage = maxDeltaDamage
	else if(deltaDamage < -maxDeltaDamage) then
		deltaDamage = -maxDeltaDamage
	end if
	
	DamageNew = DamageOld + deltaDamage
		
	if(DamageNew>MaxDamage )  DamageNew = maxDamage
	if(DamageNew<0.0)  DamageNew = 0.0
	
	DamageNew = dampFac*DamageNew + (1.0d0 - dampFac)*DamageOld
	
	
!~ 	write(0,*) "DamageNew = " , DamageNew   

	AlphaNew = 0.0d0 !! just for compatibility with ifort
       
end subroutine

subroutine GetSensibility(Wsens,XLL,Sol1,nDim,NodElT,iDofT,iShiftW)

    use funcAux
    use materialsModelsLib
    use sharedFunc

    integer:: NGP,iSimplex,iBu,pOrder, iDofVec
    parameter (NGP=1,iSimplex=0,iBu=0,pOrder=1,iDofVec=3)

    integer , intent(in) :: nDim,NodElT,iDofT,iShiftW
	real*8 , intent(in) :: XLL(*),Sol1(*)
	real*8 , intent(out) :: Wsens(Ndim)
    Real*8 :: det ,  W(NGP),Psi(nDim,NGP),Jac(nDim,nDim), GradW(nDim,nDim)
    Real*8 :: dPhi_L(nDim, NodElt), dPhi_G(nDim,NodElT),Phi(NodElT)
    Real*8 :: detF,C(nDim,nDim),F(nDim,nDim), FT(nDim,nDim), FinvT(nDim,nDim)
    Real*8 :: SolW(NodElt*iDofVec)
	integer :: i , ip, ipShiftW ,j
	Real*8 :: Udummy(ndim)
	
	W = 0.0d0
	F = 0.0d0
	Call GaussRule (Psi,W,nDim,NGP,iSimplex)

    do i=1,NodElT
        ipShiftW = (i-1) * iDofT + iShiftW
        ip = (i-1) * iDofVec
        do j=1,Ndim		
            SolW(ip + j)  = Sol1 (ipShiftW + j )
        enddo
    enddo
    
    Call LocalShapeDer(NodELT,Ndim,dPhi_L,PSI(1,NGP),pOrder,iBu) 
    Call Jacobian(Jac,det,Ndim,NodElT,XLL,dPhi_L)
    Call GlobalShapeDer(NodELT,Ndim,dPhi_L,dPhi_G,Jac)
    Call ShapeF (NodELT,Ndim,Phi,PSI(1,NGP),pOrder,iBu)

	call calcFC(Wsens,GradW,F,FT,FinvT,C,detF,SolW,Phi,dPhi_G,Ndim,NodElT,iDofVec) !! just for GradW , all the rest is dummy
   
end subroutine


subroutine IntegralSensibility(Wsens,Integral,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftW)

    use funcAux
    use materialsModelsLib
    use sharedFunc

    integer:: NGP,iSimplex,iBu,pOrder, iDofVec
    Logical :: ifAnalytical
    parameter (NGP=1,iSimplex=0,iBu=0,pOrder=1,iDofVec=3,ifAnalytical = .False.)

    integer , intent(in) :: constLaw,nDim,NodElT,iDofT,iShiftV, iShiftW
	real*8 , intent(in) :: XLL(*),Sol1(*), matPar(10)
	real*8 , intent(out) :: Integral, Wsens(Ndim)
    Real*8 :: det ,  W(NGP),Psi(nDim,NGP),Jac(nDim,nDim), GradU(nDim,nDim) , GradW(nDim,nDim), FS(Ndim,Ndim)
    Real*8 :: dPhi_L(nDim, NodElt), dPhi_G(nDim,NodElT),Phi(NodElT)
    Real*8 :: detF,C(nDim,nDim),F(nDim,nDim), FT(nDim,nDim), FinvT(nDim,nDim)
    Real*8 :: SolV(NodElt*iDofVec), SolW(NodElt*iDofVec), DDdummy(Ndim,Ndim,Ndim,Ndim) , DV
	integer :: i , ip, ipShiftV, ipShiftW ,j , nG 
	Real*8 :: Udummy(ndim) 
	
	W = 0.0d0
	F = 0.0d0
	Call GaussRule (Psi,W,nDim,NGP,iSimplex)


    do i=1,NodElT
        ipShiftV = (i-1) * iDofT + iShiftV
        ipShiftW = (i-1) * iDofT + iShiftW
        ip = (i-1) * iDofVec

        do j=1,Ndim		
            SolV(ip + j)  = Sol1 (ipShiftV + j )
            SolW(ip + j)  = Sol1 (ipShiftW + j )
        enddo
    enddo
    
!~     call printVec(SolV)
!~     call printVec(Sol1(1:NodElt*iDofT))
    
    Integral = 0.0d0
    do nG = 1 , NGP
    
    Call LocalShapeDer(NodELT,Ndim,dPhi_L,PSI(1,nG),pOrder,iBu) 
    Call Jacobian(Jac,det,Ndim,NodElT,XLL,dPhi_L)
    Call GlobalShapeDer(NodELT,Ndim,dPhi_L,dPhi_G,Jac)
    Call ShapeF (NodELT,Ndim,Phi,PSI(1,nG),pOrder,iBu)

    DV=det*W(nG)
	
!~ 	call calcFC(Wsens,GradW,F,FT,FinvT,C,detF,SolW,Phi,dPhi_G,Ndim,NodElT,iDofVec) !! just for GradW , all the rest is dummy
    F = 0.0d0
	call calcFC(Udummy,GradU,F,FT,FinvT,C,detF,SolV,Phi,dPhi_G,Ndim,NodElT,iDofVec) 
    call calcFSDD(DDdummy,FS,C,F,Ndim,matPar,constLaw,ifAnalytical,0.0d0)
		
	Integral = Integral + dot_product2(FS,GradW) * DV 
	
	end do
	
!~ 	FS = (1.0d0/detF) * (matmul(FS,FT))
!~ 	
!~ 	Integral = (FS(1,1)-FS(2,2))**2.0d0 + (FS(1,1)-FS(3,3))**2.0d0 + (FS(3,3)-FS(2,2))**2.0d0
!~ 	Integral = Integral + 6.0d0*(FS(1,2)**2.0d0 + FS(3,2)**2.0d0 + FS(1,3)**2.0d0 )
!~ 	Integral = dsqrt(0.5d0*Integral)
!~ 	Integral = detF

	WSens = 0.0d0 !! just for compatibility with ifort
	 
end subroutine

subroutine SensibilityBased(damageNew,AlphaNew,XLL,Sol1,EvolPar,matPar,damageOld,alphaOld,nDim,&
							constLaw,NodElT,iDofT,iShiftV,iShiftW)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftW ! all integers
    Real*8, intent(in) :: XLL(*),Sol1(*),EvolPar(:),matPar(:), damageOld, alphaOld
    Real*8, intent(out) :: damageNew,AlphaNew

    Real*8 :: F, vardamage, deltaDamage, dampFac, ScaleFac , W(Ndim)
    Real*8 , parameter :: MaxDamage = 0.98d0 ,  maxDeltaDamage = 0.2d0
    integer :: opSens
    Real*8 ,save :: accFac = 0.5d0 
    
    dampFac = EvolPar(1)
    ScaleFac = EvolPar(2)
    opSens = nint(EvolPar(3))
    
    select case(opSens)
		case(1) 
			call GetSensibility(W,XLL,Sol1,nDim,NodElT,iDofT,iShiftW)    
			deltaDamage = dsign(sqrt(dot_product(W,W)),W(3))/ScaleFac
		case(0)
			call IntegralSensibility(W,deltaDamage,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV+6,iShiftW)
			AlphaNew = deltaDamage
			deltaDamage = deltaDamage/ScaleFac 
    end select
    
	write(0,*) "DeltaDamageNew = " , deltaDamage
	
	if(deltaDamage > maxDeltaDamage) then 
		deltaDamage = maxDeltaDamage
	else if(deltaDamage < -maxDeltaDamage) then
		deltaDamage = -maxDeltaDamage
	end if
	
	DamageNew = DamageOld + deltaDamage
	
!~ 	accFac = 0.5d0*(accFac + dabs(deltaDamage))
!~ 	write(0,*) "accFac = " , accFac
	
	if(DamageNew>MaxDamage )  DamageNew = maxDamage
	if(DamageNew<0.0)  DamageNew = 0.0
	
	DamageNew = dampFac*DamageNew + (1.0d0 - dampFac)*DamageOld
	
	
	write(0,*) "DamageNew = " , DamageNew   
    
end subroutine

!!! No exponential, just increment as a derivative
subroutine MinimizationFunc3(damageNew,alphaNew,XLL,Sol0,Sol1,EvolPar,matPar, &
									damageOld,alphaOld,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU,iShiftVp,iShiftVm)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU,iShiftVp,iShiftVm ! all integers
    Real*8, intent(in) :: XLL(*),Sol0(*),Sol1(*),EvolPar(:),matPar(:), damageOld, alphaOld
    Real*8, intent(out) :: damageNew,alphaNew

    Real*8 :: F, Fp, Fm , weight(4) , dF, vardamage, deltaDamage, dampFac, ScaleFac
    Real*8 , parameter :: MaxDamage = 0.99d0 ,  maxDeltaDamage = 0.4d0
        
    weight(1:4) = EvolPar(1:4)
    vardamage = EvolPar(5)
    dampFac = EvolPar(6)
    ScaleFac = EvolPar(7)
   
    call CostFunc(F,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
    call CostFunc2(Fp,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftVp,iShiftU)                
    call CostFunc2(Fm,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftVm,iShiftU)
    
    dF = 0.5d0*(Fp - Fm)/vardamage
!~ 	dF = (Fp - F)/vardamage
	
	
	deltaDamage = -dF/ScaleFac
	write(0,*) "DeltaDamageNew = " , deltaDamage
	
	if(deltaDamage > maxDeltaDamage) then 
		deltaDamage = maxDeltaDamage
	else if(deltaDamage < -maxDeltaDamage) then
		deltaDamage = -maxDeltaDamage
	end if
	
	DamageNew = DamageOld + deltaDamage
	
	if(DamageNew>MaxDamage)  DamageNew = maxDamage
	if(DamageNew<0.0)  DamageNew = 0.0
	
	DamageNew = dampFac*DamageNew + (1.0d0 - dampFac)*DamageOld
	write(0,*) "DamageNew = " , DamageNew
	
	AlphaNew = 0.0d0 !! just for compatibility with ifort
		
end subroutine

!! Increment of Alpha
subroutine MinimizationFunc4(damageNew,alphaNew,XLL,Sol0,Sol1,EvolPar,matPar, &
									damageOld,alphaOld,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU,iShiftVp,iShiftVm)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU,iShiftVp,iShiftVm ! all integers
    Real*8, intent(in) :: XLL(*),Sol0(*),Sol1(*),EvolPar(:),matPar(:), damageOld, alphaOld
    Real*8, intent(out) :: damageNew,alphaNew

    Real*8 :: F, Fp, Fm , weight(4) , dF, vardamage, dampFac, ScaleFac
    Real*8 , parameter :: dinf = 0.99
        
    weight(1:4) = EvolPar(1:4)
    vardamage = EvolPar(5)
    dampFac = EvolPar(6)
	ScaleFac = EvolPar(7)
    
!~     call CostFunc(F,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
    call CostFunc(Fp,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftVp,iShiftU)                
    call CostFunc(Fm,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftVm,iShiftU)
    
    dF = 0.5d0*(Fp - Fm)/vardamage
	
	alphaNew = alphaOld - (dF/ScaleFac)
	if(alphaNew<0.0d0) alphaNew = 0.0d0
	write(0,*) "alphaOld = , alphaNew = " , alphaOld , alphaNew
	
	DamageNew = dinf*(1.0d0 - dexp(-alphaNew))
	
	DamageNew = dampFac*DamageNew + (1.0d0 - dampFac)*DamageOld
	
	write(0,*) "DamageNew exp = " , DamageNew
	
	AlphaNew = 0.0d0 !! just for compatibility with ifort
		
end subroutine

subroutine derCostFunc(Dcost,XLL,Sol0,Sol1,weight,matPar,varDamage,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU,iShiftVp)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU,iShiftVp ! all integers
    Real*8, intent(in) :: XLL(*),Sol0(*),Sol1(*),matPar(:),weight(:), varDamage
    Real*8, intent(out) :: Dcost

    Real*8 :: energyV, energyU , U(ndim), V(ndim), energyVp, Vp(ndim), dV(ndim), dEnergyV, DcostContrib(2)
        
    call calcEnergyAndDisp(V,energyV,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
	call calcEnergyAndDisp(U,energyU,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftU)
	call calcEnergyAndDisp(Vp,energyVp,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftVp)                
	
	dV = (1.0d0/varDamage) * (Vp - V)
	dEnergyV = (1.0d0/varDamage) * (energyVp - energyV)
	
	DcostContrib(1) = dot_product(V-U,dV)
	DcostContrib(2) = (energyV-energyU)*dEnergyV
	
	Dcost = 2.0*dot_product(DcostContrib,weight(1:2))
	
!~     write(6,*) "CostContrib = " ,  costContrib
!~     write(6,*) "Total Cost = " ,  cost   
      
end subroutine

subroutine CostFunc(cost,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU ! all integers
    Real*8, intent(in) :: XLL(*),Sol0(*),Sol1(*),matPar(:),weight(:)
    Real*8, intent(out) :: cost

    Real*8 :: energyV, energyU , U(ndim), V(ndim), energyV0, V0(ndim), costContrib(4)
    
    call calcEnergyAndDisp(V,energyV,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
	call calcEnergyAndDisp(U,energyU,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftU)
	call calcEnergyAndDisp(V0,energyV0,XLL,Sol0,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
	
	costContrib(1) = dot_product(V-U,V-U)
	costContrib(2) = (energyV-energyU)**2.0d0
	costContrib(3) = dot_product(V-V0,V-V0)
	costContrib(4) = (energyV-energyV0)**2.0d0

!~ 	costContrib(3) = 0.0d0
!~ 	costContrib(4) = 0.0d0
	
	cost = dot_product(costContrib,weight(1:4))
	
!~     write(0,*) "CostContrib = " ,  costContrib(2)/costContrib(1)
!~     write(0,*) "Total Cost = " ,  cost   
      
end subroutine

subroutine CostFunc2(cost,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU ! all integers
    Real*8, intent(in) :: XLL(*),Sol0(*),Sol1(*),matPar(:),weight(:)
    Real*8, intent(out) :: cost

    Real*8 :: energyV, energyU , U(ndim), V(ndim), energyV0, V0(ndim), costContrib(4) , normU , normV0 
        
    Real*8 , parameter :: tol =1.0e-12
    
    call calcEnergyAndDisp(V,energyV,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
	call calcEnergyAndDisp(U,energyU,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftU)
	call calcEnergyAndDisp(V0,energyV0,XLL,Sol0,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
	
	normU = dsqrt(dot_product(U,U))
	normV0 = dsqrt(dot_product(V0,V0))
	
	costContrib = 0.0d0
	
!~ 	costContrib(1) = dlog(1.0 + dot_product(V-U,V-U)) /dlog(1.0d0 + (normU*normU))
!~ 	costContrib(2) = dlog(1.0 + (energyV-energyU)**2.0d0) / dlog(1.0d0 + energyU*energyU) 
!~ 	costContrib(3) = 0.0d0
!~ 	costContrib(4) = 0.0d0
	
	if(normU>tol) costContrib(1) = dot_product(V-U,V-U)/(normU*normU)
	if(energyU>tol) costContrib(2) = ((energyV-energyU)/energyU)**2.0d0
	if(normV0>tol) costContrib(3) = dot_product(V-V0,V-V0)/(normV0*normV0)
	if(energyV0>tol) costContrib(4) = ((energyV-energyV0)/energyV0)**2.0d0
	
!~ 	costContrib(1) = dot_product(V-U,V-U)/(normU*normU)
!~ 	costContrib(2) = ((energyV-energyU)/energyU)**2.0d0
!~ 	costContrib(3) = dot_product(V-V0,V-V0)/(normV0*normV0)
!~ 	costContrib(4) = ((energyV-energyV0)/energyV0)**2.0d0
	
	cost = dot_product(costContrib,weight(1:4))
	
!~     write(6,*) "CostContrib = " ,  costContrib
!~     write(6,*) "Total Cost = " ,  cost   
      
end subroutine


subroutine MinimizationFunc(damageNew,alphaNew,XLL,Sol0,Sol1,EvolPar,matPar, &
									damageOld,alphaOld,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU,iShiftVp,iShiftVm)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU,iShiftVp,iShiftVm ! all integers
    Real*8, intent(in) :: XLL(*),Sol0(*),Sol1(*),EvolPar(:),matPar(:), damageOld, alphaOld
    Real*8, intent(out) :: damageNew,alphaNew

    Real*8 :: F, Fp, Fm , weight(4) , dF, ddF, vardamage, deltaDamage, dampFac, dFsqrt
    Real*8 , parameter :: MaxDamage = 0.99d0 ,  maxDeltaDamage = 0.2d0 , tol = 1.0d0-4
        
    weight(1:4) = EvolPar(1:4)
    vardamage = EvolPar(5)
    dampFac = EvolPar(6)
   
    call CostFunc(F,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
    call CostFunc(Fp,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftVp,iShiftU)                
    call CostFunc(Fm,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftVm,iShiftU)
    
    ddF = (Fp + Fm - 2.0d0*F)/(vardamage*vardamage)
    dF = 0.5d0*(Fp - Fm)/vardamage
	
	deltaDamage = 0.0d0
	if(dabs(F)>tol) deltaDamage = -F/dF

	deltaDamage = dampFac*deltaDamage
	
	
	if(deltaDamage > maxDeltaDamage) then 
		deltaDamage = maxDeltaDamage
	else if(deltaDamage < -maxDeltaDamage) then
		deltaDamage = -maxDeltaDamage
	end if
	
	DamageNew = DamageOld + deltaDamage
	
	if(DamageNew>MaxDamage)  DamageNew = maxDamage
	if(DamageNew<0.0)  DamageNew = 0.0
	
	AlphaNew = 0.0d0 !! just for compatibility with ifort
		
end subroutine


subroutine MinimizationFunc2(damageNew,alphaNew,XLL,Sol0,Sol1,EvolPar,matPar, &
									damageOld,alphaOld,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU,iShiftVp,iShiftVm)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU,iShiftVp,iShiftVm ! all integers
    Real*8, intent(in) :: XLL(*),Sol0(*),Sol1(*),EvolPar(:),matPar(:), damageOld, alphaOld
    Real*8, intent(out) :: damageNew,alphaNew

    Real*8 :: F, Fp, Fm , weight(4) , dF, ddF, vardamage, deltaDamage, dampFac, dFsqrt , damageNew1,damageNew2
    Real*8 , parameter :: MaxDamage = 0.99d0 ,  maxDeltaDamage = 0.8d0 , tol = 0.0e-8
        
    weight(1:4) = EvolPar(1:4)
    vardamage = EvolPar(5)
    dampFac = EvolPar(6)
   
    call CostFunc2(F,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
    call CostFunc2(Fp,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftVp,iShiftU)                
    call CostFunc2(Fm,XLL,Sol0,Sol1,weight,matPar,nDim,constLaw,NodElT,iDofT,iShiftVm,iShiftU)
    
    ddF = (Fp + Fm - 2.0d0*F)/(vardamage*vardamage)
    dF = 0.5d0*(Fp - Fm)/vardamage

!~     call derCostFunc(dF,XLL,Sol0,Sol1,weight,matPar,varDamage,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU,iShiftVp ) 
	
!~ 	deltaDamage = 0.0d0
!~ 	if(dabs(F)>tol) deltaDamage = -F/dF
	
!~ 	deltaDamage = -dF/ddF
!~ 	deltaDamage = -2.0d0*F*dF/(2.0d0*dF*dF - F*ddF) !! Halley's Method
!~ 	deltaDamage = deltaDamage - 0.5d0*F*F*ddF/(dF**3.0d0) !! Chebyschev's Method
	
!~ 	deltaDamage = dampFac*deltaDamage
	
!~ 	write(6,*) "deltaDamage = " , deltaDamage
!~ 	write(6,*)  "F=" , F, "Fp=" , Fp, "Fm=" , Fm,  "dF = ", dF, "ddF = ", ddF , "DeltaDamage = " , deltaDamage
		
!~ 	if(dabs(F*ddF)<dF*dF) then
!~ 		deltaDamage = -F/dF
!~ 		write(0,*) "choosen NR"
!~ 	else 
!~ 		deltaDamage = -dF/ddF
!~ 		write(0,*) "choosen Newton"
!~ 	end if
		
	
	deltaDamage = -F/dF
	
	if(deltaDamage > maxDeltaDamage) then 
		deltaDamage = maxDeltaDamage
	else if(deltaDamage < -maxDeltaDamage) then
		deltaDamage = -maxDeltaDamage
	end if

	DamageNew1 = DamageOld + deltaDamage
	
	if(DamageNew1>MaxDamage)  DamageNew1 = maxDamage
	if(DamageNew1<0.0)  DamageNew1 = 0.0

	DamageNew1 = DamageNew1*dampFac + DamageOld*(1.0d0 - dampFac)
	
	deltaDamage = -dF/ddF
	
	if(deltaDamage > maxDeltaDamage) then 
		deltaDamage = maxDeltaDamage
	else if(deltaDamage < -maxDeltaDamage) then
		deltaDamage = -maxDeltaDamage
	end if

	DamageNew2 = DamageOld + deltaDamage
	
	if(DamageNew2>MaxDamage)  DamageNew2 = maxDamage
	if(DamageNew2<0.0)  DamageNew2 = 0.0

	DamageNew2 = DamageNew2*dampFac + DamageOld*(1.0d0 - dampFac)

	damageNew = 0.5d0*(damageNew2 + damageNew1)
!~ 	if(dabs(damageNew2 - damageOld)<tol ) then
!~ 		damageNew = damageNew1
!~ 	else if(dabs(damageNew1 - damageOld)<tol ) then
!~ 		damageNew = damageNew2
!~ 	else 
!~ 		damageNew = min(damageNew2,damageNew1)
!~ 	end if

	AlphaNew = 0.0d0 !! just for compatibility with ifort
end subroutine


subroutine calcEnergyAndDisp(U,energy,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShift)

    use funcAux
    use materialsModelsLib
    use sharedFunc

    integer:: NGP,iSimplex,iBu,pOrder, iDofVec
    parameter (NGP=1,iSimplex=0,iBu=0,pOrder=1,iDofVec=3)

    integer , intent(in) :: constLaw,nDim,NodElT,iDofT,iShift
	real*8 , intent(in) :: XLL(*),Sol1(*), matPar(10)
	real*8 , intent(out) :: energy, U(ndim)
    Real*8 :: det ,  W(NGP),Psi(nDim,NGP),Jac(nDim,nDim), GradU(nDim,nDim)
    Real*8 :: dPhi_L(nDim, NodElt), dPhi_G(nDim,NodElT),Phi(NodElT)
    Real*8 :: detF,C(nDim,nDim),F(nDim,nDim), FT(nDim,nDim), FinvT(nDim,nDim)
    Real*8 :: SolShift(NodElt*iDofVec)
	integer :: i , ip, ipShift,j,jp,jpShift
	
	W = 0.0d0
	F = 0.0d0
	Call GaussRule (Psi,W,nDim,NGP,iSimplex)

    do i=1,NodElT

        ipShift = (i-1) * iDofT + iShift
        ip = (i-1) * iDofVec

        do j=1,Ndim
        
			jpShift = ipShift + j
			jp = ip + j
			
            SolShift(jp)  = Sol1 (jpShift)
        enddo
    enddo
    
    Call LocalShapeDer(NodELT,Ndim,dPhi_L,PSI(1,NGP),pOrder,iBu) 
    Call Jacobian(Jac,det,Ndim,NodElT,XLL,dPhi_L)
    Call GlobalShapeDer(NodELT,Ndim,dPhi_L,dPhi_G,Jac)
    Call ShapeF (NodELT,Ndim,Phi,PSI(1,NGP),pOrder,iBu)
	
	F = 0.0d0
    call calcFC(U,GradU,F,FT,FinvT,C,detF,SolShift,Phi,dPhi_G,Ndim,NodElT,iDofVec)
    call strainEnergy(energy,C,matPar,constLaw)

end subroutine


subroutine energyExponential(damageNew,alphaNew,XLL,Sol1,EvolPar,matPar, &
							damageOld,alphaOld,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
							
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU ! all integers
    Real*8, intent(in) :: XLL(*),Sol1(*),EvolPar(:),matPar(:),alphaOld, damageOld
    Real*8, intent(out) :: damageNew,alphaNew

    Real*8 :: energyV, energyU , dinf, beta
    Real*8 :: vDummy(ndim)
    
    dinf = EvolPar(1)
    beta = EvolPar(2)
    
!~     write(0,*) "calculating energyV"
    call calcEnergyAndDisp(vDummy,energyV,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
!~ 	 write(0,*) "calculating energyU"
	call calcEnergyAndDisp(vDummy,energyU,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftU)
        
    alphaNew = alphaOld + energyU - energyV
    if(alphaNew<0.0d0) alphaNew = 0.0d0
   
    damageNew = dinf*(1.0d0 - dexp(-alphaNew/beta)) 
    
!~ 	write(6,*) "DeltaDamage=",damageNew-damageOld
!~     write(6,*) "Damage=",damageNew,"alpha=",alphaNew
!~     write(6,*) "EnergyU=",energyU,"energyV=",energyV
!~     
    
end subroutine

subroutine energyExponentialNorm(damageNew,alphaNew,XLL,Sol0,Sol1,EvolPar,matPar, &
									damageOld,alphaOld,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU ! all integers
    Real*8, intent(in) :: XLL(*),Sol0(*),Sol1(*),EvolPar(:),matPar(:), damageOld
    Real*8, intent(out) :: damageNew,alphaNew

    Real*8 :: energyV, energyU , energyV0, r , beta,dinf, dinf0, dinfF, beta0, betaF, alphaOld
    Real*8 :: vDummy(ndim)
    
    dinf0 = EvolPar(1)
    dinfF = EvolPar(2)
    beta0 = EvolPar(3)
    betaF = EvolPar(4)
    
    dinf = ((dinfF-damageOld)*dinf0 + damageOld*dinfF)/dinfF 
    
    beta = ((dinf-damageOld)*beta0 + damageOld*betaF)/dinf 
    
    if(damageOld == 0.0d0) then
		alphaOld = 0.0d0
    else 
		alphaOld = -beta*dlog(1.0 - damageOld/dinf)
    end if 
    
    call calcEnergyAndDisp(vDummy,energyV,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
!~ 	call calcEnergy(energyV0,XLL,Sol0,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
	call calcEnergyAndDisp(vDummy,energyU,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftU)
    
!~     alphaNew = alphaOld + (energyU - energyV) + dabs(energyV-energyV0)
	if( (energyU - energyV) > 0.0d0) then
		alphaNew = alphaOld + (energyU - energyV)
    else 
		alphaNew = alphaOld + 10.0d0*(energyU - energyV)
	end if 
	
    if(alphaNew<0.0d0) alphaNew = 0.0d0
   
!~     damageNew = dinf*(1.0d0 - dexp(-alphaNew/beta))**(1.0d0 - damageNew)
	 
    damageNew = dinf*(1.0d0 - dexp(-alphaNew/beta))
    
    if(damageNew > 0.95 ) write(6,*) "damageNew = " , damageNew
      
end subroutine

subroutine energyLinear(damageNew,alphaNew,XLL,Sol1,EvolPar,matPar,damageOld,alphaOld,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU ! all integers
    Real*8, intent(in) :: XLL(*),Sol1(*),EvolPar(:),matPar(:),alphaOld , damageOld 
    Real*8, intent(out) :: damageNew,alphaNew 

    Real*8 :: energyV, energyU , kappa, MaxDeltaDamage
    Real*8 :: vDummy(ndim)
    
    
    kappa = EvolPar(1)
    MaxDeltaDamage = EvolPar(2)
    
    call calcEnergyAndDisp(vDummy,energyV,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
	call calcEnergyAndDisp(vDummy,energyU,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftU)
   
    alphaNew = alphaOld + energyU - energyV
    damageNew = kappa*alphaNew
    
	if(damageNew > 0.99d0 ) then
		damageNew = 0.99d0
	else if(damageNew < 0.0d0 ) then
		damageNew= 0.0d0
	else if(damageNew - damageOld > MaxDeltaDamage) then 
		damageNew = damageOld + MaxDeltaDamage
	else if(damageNew - damageOld < -MaxDeltaDamage) then 
		damageNew = damageOld - MaxDeltaDamage
	end if

	write(6,*) "DeltaDamage=",damageNew-damageOld
    write(6,*) "Damage=",damageNew,"alpha=",alphaNew
    write(6,*) "EnergyU=",energyU,"energyV=",energyV
    
end subroutine

subroutine energyNormal(damageNew,alphaNew,XLL,Sol1,EvolPar,matPar, &
							damageOld,alphaOld,nDim,constLaw,NodElT,iDofT,iShiftV,iShiftU)
							
	integer, intent(in) :: Ndim,iDofT,NodELT,constLaw,iShiftV,iShiftU ! all integers
    Real*8, intent(in) :: XLL(*),Sol1(*),EvolPar(:),matPar(:), damageOld
    Real*8, intent(out) :: damageNew,alphaNew

    Real*8 :: energyV, energyU , dinf, beta, dm,alphaOld, mu, fac
    Real*8 :: vDummy(ndim)
    
    dinf = EvolPar(1)
    dm = EvolPar(2)
    mu = EvolPar(3)
    
!~     mu = 2.0d0*beta*dsqrt(-dlog(dm/dinf))
    beta = 0.5d0*mu/dsqrt(-dlog(dm/dinf))
    alphaOld = mu + beta*dsqrt(-dlog(damageOld/dinf))
    
!~     write(0,*) "calculating energyV"
    call calcEnergyAndDisp(vDummy, energyV,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftV)                
!~ 	 write(0,*) "calculating energyU"
	call calcEnergyAndDisp(vDummy, energyU,XLL,Sol1,matPar,nDim,constLaw,NodElT,iDofT,iShiftU)
        
    alphaNew = alphaOld + energyU - energyV
!~     write(0,*) alphaNew
!~     if(alphaNew>mu) alphaNew = mu
   
	fac = (alphaNew - mu)/beta
    damageNew = dinf*dexp(-fac*fac) 
    
!~ 	write(6,*) "DeltaDamage=",damageNew-damageOld
!~     write(6,*) "Damage=",damageNew,"alpha=",alphaNew
!~     write(6,*) "EnergyU=",energyU,"energyV=",energyV
!~     
    
end subroutine


end module
