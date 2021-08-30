!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity 
! Still in development
!     ------------------------------------------------------------------

!~ 
Subroutine PostProcSensib(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, &
				Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use materialsModelsLib
	use sharedFunc
	use damageCriteria
	
	implicit none
	
	! for compatibility is Real*8, integer, integer. Det to me is local, not shared in common
	real*8 :: dummy
	integer:: LengthParam,LengthJParam
	COMMON dummy,LengthParam,LengthJParam

	!	PARAMETERS
	
	integer :: iSimplex, NGPMAX, iBu , NParamField , OElemOutput
	parameter (iBu = 0, NGPMAX = 4, iSimplex = 0, NParamField=2, OElemOutput = 500) ! No Bubbles and Linear, ! iSimplex = 0 -> The element is a Simplex 
!~ 

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

	!   =====   END ARGUMENTS  =======
         
	Integer :: nG,i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp, ip , ipp! counters 
	Integer :: NDimE,NodG,NodEL, NparamElem
	Parameter( NdimE = 3) ! this code only makes sense in 3D  
	Real*8 :: Det
	Real*8 :: Psi(NdimE,NGPMAX), W(NGPMAX), Phi(NodElT), & 
			  dPhi_G(NdimE,NodElT), dPhi_L(NdimE,NodElT), Jac(NDimE,NDimE)
	Real*8 :: U(NdimE), GradU(NdimE,NdimE) , F(NdimE,NdimE), FT(NdimE,NdimE) , FinvT(NdimE,NdimE) , C(NdimE,NdimE) , detF
	Real*8 :: FS(NdimE,NdimE), DD(NdimE,NdimE,NdimE,NdimE)
	Real*8 :: MatPar(10) , Sol1V(NdimE*NodElT) , Sol1W(NdimE*NodElT) , Sol1U(NdimE*NodElT), Sol1VU(NdimE*NodElT),VField(4)
	Real*8 :: normL2_ref , normEnergy_ref , weight1 , weight2 , GradW(NdimE,NdimE) , tempVolModulus
	integer :: iShiftV , iShiftW , iShiftU , ipShiftV , ipShiftW , ipShiftU 
	integer :: NGP,pOrder , IisPeriodic ,iShiftParam, constLaw , NFields, CodeField(5) , PosField(5) !! maximum 5 at now
	logical :: isPeriodic, ifAnalytical

	matPar = 0.0d0
	
	NparamElem = nint(CommonPar(1))
	pOrder = nint(commonPar(2))
	ifAnalytical = .false.
	if(nint(commonPar(3)) == 1) ifAnalytical =.true.
	constLaw = nint(CommonPar(4)) 	
	matPar(1:3) = CommonPar(5:7)
	weight1 = CommonPar(8)
	weight2 = CommonPar(9)
	normL2_ref = CommonPar(10)
	normEnergy_ref = CommonPar(11)
	NFields = nint(CommonPar(12))
	do i = 1, NFields
		CodeField(i) = nint(CommonPar( 12 + (i-1)*NParamField + 1) )
		PosField(i)  = nint(CommonPar( 12 + (i-1)*NParamField + 2) )
	end do
	
	VField = 0.0d0
	NGP = 1
	nG = 1
	NodG = 4
	
	weight1 = weight1/normL2_ref
	weight2 = weight2/normEnergy_ref
		
!~ 	write(0,*) "w1 = " , weight1 , " , w2 = ", weight2
		
	W = 0.0d0		
	Call GaussRule (Psi,W,NdimE,NGP,iSimplex)

	AE = 0.0d0
	
	Call LocalShapeDer(NodG,NdimE,dPhi_L,PSI(1,nG),pOrder,iBu) ! PSI(1,nG) is a pointer to the nG^th column 
	Call Jacobian(Jac,Det,NdimE,NodG,XLL,dPhi_L)
	Call GlobalShapeDer(NodG,NdimE,dPhi_L,dPhi_G,Jac)
	Call ShapeF (NodELT,NdimE,Phi,PSI(1,nG),pOrder,iBu)
	
	iShiftV = 0
	iShiftW = 3
	iShiftU = 6
	do i=1,NodElT
        ipShiftV = (i-1) * iDofT + iShiftV
        ipShiftW = (i-1) * iDofT + iShiftW
		ipShiftU = (i-1) * iDofT + iShiftU
        ip = (i-1) * NdimE

        do j=1,Ndim		
            Sol1V(ip + j)  = Sol1 (ipShiftV + j )
            Sol1W(ip + j)  = Sol1 (ipShiftW + j )
            Sol1U(ip + j)  = Sol1 (ipShiftU + j )
        enddo
    enddo
	
!~ 	Sol1VU = Sol1V - Sol1U !! u - ud formulation
	Sol1VU = Sol1U - Sol1V !! u - ud formulation
!~ 	write(6,*) "Sol1W"
!~ 	call printVec(Sol1W)
	
	
!~!! 	call setDamage2Matpar(matPar,Param,constLaw,NparamElem,nG)
	matPar(4) =0.0d0 !! Just for safety
	if(constLaw == 9) constLaw = 3
	if(constLaw == 10) constLaw = 4
	
	do i = 1 , NFields
		select case(CodeField(i))
			case(1) !! volume
				VField(i) = det*W(nG)
			case(2) !! damage
				VField(i) = Param(LengthParam + (nG-1)*NParamElem + 1 )
			case(3) !! ||V - U||^2 
				F = 0.0d0
				call calcFC(U,GradU,F,FT,FinvT,C,detF,Sol1VU,Phi,dPhi_G,NdimE,NodG,Ndim)
				VField(i) = dot_product(Sol1VU,Sol1VU)
				J_km = J_k
!~ 				write(0,*) VField(i) , "case 3"
			case(4) !! energy V - U
				F = 0.0d0
				call calcFC(U,GradU,F,FT,FinvT,C,detF,Sol1VU,Phi,dPhi_G,NdimE,NodG,Ndim)
				call calcFSDD(DD,FS,C,F,NdimE,matPar,constLaw,ifAnalytical,0.0d0)
				call strainEnergy(VField(i),C,matPar,constLaw)
!~ 				write(0,*) VField(i) , "case 4"
			case(5) !! Functional = w1*||V-U||^2/normL2_ref + w2*energy(V-U)/normEnergy_ref
				F = 0.0d0 
				call calcFC(U,GradU,F,FT,FinvT,C,detF,Sol1VU,Phi,dPhi_G,NdimE,NodG,Ndim)
				call calcFSDD(DD,FS,C,F,NdimE,matPar,constLaw,ifAnalytical,0.0d0)
				call strainEnergy(VField(i),C,matPar,constLaw)
				VField(i) = weight2*VField(i)
				VField(i) = VField(i) + weight1*dot_product(Sol1VU,Sol1VU)
				J_k = J_k + det*W(nG)*VField(i)
				
				Idam = Idam + det*W(nG)*Param(1)
				measOmega = measOmega + det*W(nG)
				
			case(6) !! Integral
				F = 0.0d0
				call calcFC(U,GradU,F,FT,FinvT,C,detF,Sol1V,Phi,dPhi_G,NdimE,NodG,Ndim) 
				tempVolModulus = matPar(3)
				MatPar(3) = 0.0d0
				call calcFSDD(DD,FS,C,F,NdimE,matPar,constLaw,ifAnalytical,0.0d0) 
				matPar(3) = tempVolModulus
				call calcFC(U,GradW,F,FT,FinvT,C,detF,Sol1W,Phi,dPhi_G,NdimE,NodG,Ndim) !! Just for GradW
				VField(i) = dot_product2(FS,GradW)*det*W(nG)
				if(dabs(VField(i))>Ikmax) Ikmax = dabs(VField(i)) 
			case default
				write(0,*) "Code Field not found"
				VField(i) = 0
		end select
		
		Param(LengthParam + (nG-1)*NParamElem + PosField(i) ) = VField(i)
	end do
	
!~ 	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
		
End Subroutine



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine PostProcSensibS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IMPLICIT NONE
	Integer :: MaxLRows,Ndim,iDofT, iAdd, ipRow , i
	Integer   Coupling(MaxLRows,MaxLRows)
	Real*8 CommonPar(*)

	Coupling=0 
	
	Coupling(1,1) = 1
	
	return
end Subroutine



