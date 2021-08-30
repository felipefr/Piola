!Last Modification : 11/03/2014
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Element for Hyperelasticity
! Still in development
!     ------------------------------------------------------------------
Subroutine AdjProbDamageSens(AE, BE, MaxLRows, XLL, NDim, iDofT, NodELT, & 
							  Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
    !     ------------------------------------------------------------------
    use funcAux
    use sharedFunc
	use damageCriteria
	
    implicit none

    !   PARAMETERS
    integer :: iSimplex, NGPMAX, iBu
    parameter (iBu = 0, NGPMAX = 11, iSimplex = 0) ! No Bubbles and Linear, ! iSimplex = 0 -> The element is a Simplex
    !~  parameter (iBu = 0, pOrder = 2,NGP = 11, iSimplex = 0) ! No Bubbles and Linear, ! iSimplex = 0 -> The element is a Simplex

    !   ===== SUBROUTINE ARGUMENTS  =======
    integer :: MaxLRows,Ndim,iDofT,NodELT ! all integers
    Real*8 :: DelT, DTm,Time ! all reals
    Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

    !   =====   END ARGUMENTS  =======

    Integer :: nG, i,j,k,l,m,n,p,q, A, B, ApRow, BpCol , App,Bqp , ApDim ! counters
    Integer :: NDimE,NodG,NodEL, NparamDamage, constLaw ! if NparamDamage = 0 , no damage
    logical :: ifAnalytical
    Parameter( NdimE = 3) ! this code only makes sense in 3D
    Real*8 :: PhiA,PhiB, dv , Det , MatPar(10)
    Real*8 :: Psi(NdimE,NGPMAX), W(NGPMAX), Phi(NodElT), dPhiA(NdimE), dPhiB(NdimE),&
        dPhi_G(NdimE,NodElT), dPhi_L(NdimE,NodElT), Jac(NDimE,NDimE) , Sol1U(NdimE*NodElt) , Sol1V(NdimE*NodElt)
        
    Real*8 :: U(NdimE), GradU(NdimE,NdimE) , F(NdimE,NdimE), FT(NdimE,NdimE) , FinvT(NdimE,NdimE) , C(NdimE,NdimE) , detF
    Real*8 :: S(NdimE,NdimE), FinvTGradU, DD(NdimE,NdimE,NdimE,NdimE) , DDGradU(NdimE,NdimE), FS(NdimE,NdimE) , & 
			 weight1 ,weight2 ,  normL2_ref , normEnergy_ref , Sol1VU(NdimE*NodElt) , Pvu(NdimE,NdimE)
    integer :: pOrder , NGP, iShiftV ,iShiftW, iShiftU , constLawUndamaged

    NodG  = NodELT
    NodEL = NodELT
    
    constLawUndamaged = 0

    ! Constitutive Parameters (commonPar(2),commonPar(3)) physical interpretation depends on the constitutive law chosen
    !~  call chooseLaw(nint(CommonPar(1)))
    constLaw = nint(CommonPar(1))
    if(constLaw == 9) constLawUndamaged = 3
    if(constLaw == 10) constLawUndamaged = 4
    matPar(1:3) = commonPar(2:4)
    ifAnalytical = .false.
    if(commonPar(5) == 1) ifAnalytical =.true.
    NparamDamage = nint(commonPar(6))
    pOrder = nint(commonPar(7))
    weight1 = CommonPar(8)
	weight2 = CommonPar(9)
	normL2_ref = CommonPar(10)
	normEnergy_ref = CommonPar(11)
	
	weight1 = weight1/normL2_ref
	weight2 = weight2/normEnergy_ref
	
	
!~ 	write(0,*) "w1 = " , weight1 , " , w2 = ", weight2
	
	iShiftV = 0
	iShiftW = 3
	iShiftU = 6 
	
	do j = 1 , NdimE
		Sol1U(j:(NodElT-1)*NdimE + j:NdimE) = Sol1(j+iShiftU:(NodElT-1)*idofT + j + iShiftU: idofT)
		Sol1V(j:(NodElT-1)*NdimE + j:NdimE) = Sol1(j+iShiftV:(NodElT-1)*idofT + j + iShiftV: idofT) 
	end do 
	
!~ 	Sol1VU = Sol1V - Sol1U !! V - U !! u - ud formulation
	Sol1VU = Sol1U - Sol1V !! U - V  !! ud - u formulation

    select case(pOrder)
        case(1)
            NGP = 4
        case(2)
            NGP = 11
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
        Call LocalShapeDer(NodELT,NdimE,dPhi_L,PSI(1,nG),pOrder,iBu) ! PSI(1,nG) is a pointer to the nG^th column
        if ( nG .eq. 1 ) then
            Call Jacobian(Jac,Det,NdimE,NodG,XLL,dPhi_L)
        endif
        Call GlobalShapeDer(NodELT,NdimE,dPhi_L,dPhi_G,Jac)
        Call ShapeF (NodELT,NdimE,Phi,PSI(1,nG),pOrder,iBu)

        call setDamage2Matpar(matPar,Param,constLaw,NParamDamage,nG) 
        		
		F = 0.0d0
        call calcFC(U,GradU,F,FT,FinvT,C,detF,Sol1VU,Phi,dPhi_G,NdimE,NodElT,NdimE)
        call calcFSDD(DD,Pvu,C,F,NdimE,matPar,constLawUndamaged,ifAnalytical,0.0d0)
		
		F = 0.0d0
        call calcFC(U,GradU,F,FT,FinvT,C,detF,Sol1V,Phi,dPhi_G,NdimE,NodElT,NdimE)
        call calcFSDD(DD,FS,C,F,NdimE,matPar,constLaw,ifAnalytical,0.0d0)

        DV=Det*W(nG)
 
        FinvTGradU = dot_product2(FinvT,GradU)
        call T4xT2(DDGradU,DD,GradU)

        Do A=1,NodEL !LoopRow
            ApRow  = (A-1)*iDofT + iShiftW
            ApDim  = (A-1)*NdimE
           
            PhiA = Phi(A)
            dPhiA=dPhi_G(:,A)

            Do p = 1 , NdimE
                App = ApRow+p

				!! u - ud formulation
				!BE(App) = BE(App) - 2.0d0 * weight1 * Sol1VU(ApDim + p) * PhiA * DV
				!BE(App) = BE(App) - weight2 * dot_product(Pvu(p,:),dPhiA) * DV
				
				BE(App) = BE(App) + 2.0d0 * weight1 * Sol1VU(ApDim + p) * PhiA * DV
				BE(App) = BE(App) + weight2 * dot_product(Pvu(p,:),dPhiA) * DV
				
				
                do B=1,NodEL ! LoopCol ! not considering the simmetry
                    BpCol  = (B-1)*iDofT + iShiftW
                    PhiB   = Phi(B)
                    dPhiB=dPhi_G(:,B)

                    do q=1,NdimE
                        Bqp = BpCol+q

						AE(App,Bqp) = AE(App,Bqp) + DDGradPhiBdotGradPhiA(DD,dPhiA,dPhiB,q,p,NdimE)*DV

                    end do ! loop Bq
                Enddo !LoopCol
            end do ! loop Ap
        Enddo !LoopRow
    EndDo !LoopGauss

	Ikmax = 0.0d0
	J_k = 0.0d0
	Idam = 0.0d0
	measOmega = 0.0d0
	    
End subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine AdjProbDamageSensS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Integer :: A,B,q,p, ApRow, BpCol, iShiftW, NodElT

    NodElT = MaxLRows / iDofT

    Coupling = 0
!~ 	Coupling(1,1) = 1

    iShiftW = 3
!~ 
    do A=1,NodElt
        ApRow = (A-1) * iDofT + iShiftW
		
		do B = 1,NodElt
			BpCol = (B-1) * iDofT + iShiftW

			do p=1,nDim
			do q=1,nDim

				Coupling (ApRow + p, BpCol + q) = 1

			enddo
            enddo
		enddo
	enddo
	
end Subroutine




