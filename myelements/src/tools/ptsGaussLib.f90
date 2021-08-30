module ptsGaussLib
	use funcAux
	implicit none
	
	type :: PtGaussClass
	
		integer :: NGP, NodG, NdimE
		real*8 , allocatable :: dPhi_G(:,:,:), Phi(:,:) , dV(:) 
	
		Contains
				
		procedure, public :: init, calcGradU, calcU
		
	end type
	
	public GaussRuleQ, setFEMtype, setNodG
	
	contains 
	
	subroutine setNodG(FEMtype,NodG)
		integer , intent(in) :: FEMtype
		integer , intent(out) :: NodG
		
		select case(FEMtype)
			case(1) !! linear triangle
				NodG = 3
			case(2) !! linear quadrangle
				NodG = 4
			case(3) !! linear quadrangle
				NodG = 3
			case default
				write(*,*) "FEMtype = " , FEMtype,  "not supported at moment"  
				stop
		end select

	end subroutine 
	
	subroutine setFEMtype(FEMtype,NodG,pOrder,NGP,iSimplex,iBubble)
		integer , intent(in) :: FEMtype
		integer , intent(out) :: NodG,pOrder,NGP,iSimplex,iBubble
		
		
		select case(FEMtype)
			case(1) !! linear triangle
				NGP = 1
				NodG = 3
				pOrder = 1
				iSimplex = 0
				iBubble = 0
		
			case(2) !! linear quadrangle
				NGP = 1
				NodG = 4
				pOrder = 1
				iSimplex = 1
				iBubble = 0
		
			case(3) !! linear triangle with 3 gauss points
				NGP = 3
				NodG = 3
				pOrder = 1
				iSimplex = 0
				iBubble = 0
		
			case default
				write(*,*) "FEMtype = " , FEMtype,  "not supported at moment"  
				stop
		end select
		
	end subroutine 

	
	subroutine GaussRuleQ(Psi,W,NdimE,NGP, iSimplex)
		integer , intent(in) :: NdimE, NGP,  iSimplex
		real*8, intent(out) :: Psi(NdimE,NGP), W(NGP)
		integer , parameter :: NGP1Dmax = 3, Ndim1D = 1
		integer :: nGi, nGj, nGk, NGP1D, nGpp
		real*8 :: Psi1D(1,NGP1Dmax), W1D(NGP1Dmax) 
	
		NGP1D = nint(real(NGP)**(1.0d0/real(NdimE)))
		
		Call GaussRule(Psi1D(:,1:NGP1D) ,W1D(1:NGP1D),Ndim1D,NGP1D,iSimplex)
		
		select case(NdimE)
			case(1)
				Psi(1,:) = Psi1D(1,1:NGP1D)
				W = W1D(1:NGP1D)
			case(2)
				do nGi = 1,NGP1D
				do nGj = 1,NGP1D
					nGpp = (nGi - 1)*NGP1D + nGj 
					Psi(1,nGpp) = Psi1D(1,nGi)
					Psi(2,nGpp) = Psi1D(1,nGj)
					W(nGpp) = W1D(nGi) * W1D(nGj)
				end do
				end do
			
			case(3)
				do nGi = 1,NGP1D
				do nGj = 1,NGP1D
				do nGk = 1,NGP1D
					nGpp = (nGi - 1)*NGP1D*NGP1D + (nGj-1)*NGP1D + nGk 
					Psi(1,nGpp) = Psi1D(1,nGi)
					Psi(2,nGpp) = Psi1D(1,nGj)
					Psi(3,nGpp) = Psi1D(1,nGk)
					W(nGpp) = W1D(nGi) * W1D(nGj) * W1D(nGk)
				end do
				end do
				end do
		end select
	
	
	end subroutine
	
	
	subroutine init(s,XLL,NodG,NdimE,NGP, pOrder,iBu, iSimplex)
		class(PtGaussClass) :: s
		integer , intent(in) :: NodG, NGP, NdimE, iBu, iSimplex, pOrder
		real*8 , intent(in) :: XLL(NdimE*NodG)
		integer :: i
		real*8 :: Det, W(NGP), Psi(NdimE,NGP), Jac(NdimE,NdimE), dPhi_L(NdimE,NodG)
			
		s%NodG = NodG
		s%NdimE = NdimE
		s%NGP = NGP 	
		
		allocate(s%dPhi_G(NdimE,NodG,NGP),s%Phi(NodG,NGP),s%dV(NGP))
			
		if(iSimplex == 0) then
			Call GaussRule (Psi,W,NdimE,NGP,iSimplex)
		else
			Call GaussRuleQ (Psi,W,NdimE,NGP,iSimplex)
		end if
		
		do i = 1 , NGP
			if(iSimplex == 0) then
				Call LocalShapeDer(NodG,NdimE,dPhi_L,PSI(1,i),pOrder,iBu)
				Call ShapeF (NodG,NdimE,s%Phi(1,i),PSI(1,i),pOrder,iBu)
			else 
				Call LocalShapeDerQ(NodG,NdimE,dPhi_L,PSI(1,i),pOrder,iBu)
				Call ShapeFQ (NodG,NdimE,s%Phi(1,i),PSI(1,i),pOrder,iBu)			 	
			end if
			
			Call Jacobian(Jac,Det,NdimE,NodG,XLL,dPhi_L)
			Call GlobalShapeDer(NodG,NdimE,dPhi_L,s%dPhi_G(1,1,i),Jac)
			s%dV(i) = W(i)*Det
			
		end do
		
	end subroutine
	
	subroutine calcGradU(s,GradU,SolU,nG)
		class(PtGaussClass) :: s
		integer, intent(in) :: nG 
		Real*8, intent(out) :: GradU(:,:) !! have NdimE x NdimE dimension
		Real*8, intent(in) :: SolU(:) !! have NdimE*NodG dimension
		integer :: i, j , e ,ep 
		
		GradU = 0.0d0
		do e = 1 , s%NodG
			ep = (e-1)*s%NdimE
			do i = 1, s%NdimE
				do j = 1, s%NdimE
					GradU(i,j) = GradU(i,j) + SolU(ep+i)*s%dPhi_G(j,e,nG)
				end do
			end do
		end do	
		
	end subroutine

	subroutine calcU(s,U,SolU,nG)
		class(PtGaussClass) :: s
		integer, intent(in) :: nG 
		Real*8, intent(out) :: U(:) !! have NdimE dimension
		Real*8, intent(in) :: SolU(:) !! have NdimE*NodG dimension 
		integer :: i, e ,ep 
		
		U= 0.0d0
		do e = 1 , s%NodG
			ep = (e-1)*s%NdimE
			do i = 1, s%NdimE
				U(i) = U(i) + SolU(ep+i)*s%Phi(e,nG)
			end do
		end do	
		
	end subroutine
	


end module
