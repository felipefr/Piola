module efemLib
	use funcAux
	implicit none

	public posProcessBeta

	contains

	Subroutine posProcessBeta(Beta,PtG,Param,matPar,damagePar,SolU,SolDeltaU, cTs, aTs, &
							damageTs, AreaS, constLaw, iFtype,SolitaryNode, NodG, NGP, iDamageParType)
	!     ------------------------------------------------------------------
		use funcAux
		use globalVariables, only : getF, NdimE, getDamagePar, getMaterial, maxDamagePar, maxMatPar
		use finiteStrainLib
		use damageNewLib
		use ptsGaussLib2
		
		implicit none
		
		real*8, intent(inout) :: Beta(NdimE)
		type(ptGaussClass), intent(in) :: PtG
		real*8 , intent(in) :: matPar(:), damagePar(:), SolU(:), SolDeltaU(:)
		real*8 , intent(in) :: cTs, aTs, damageTs, AreaS
		integer , intent(in) :: NGP, NodG, iFtype, SolitaryNode, constLaw, iDamageParType
		real*8 :: Param(*)
		

	!~ 	!   =====   END ARGUMENTS  =======
			 
		Integer :: nG ! counters 
		Real*8 :: dV , Det
		Real*8 :: GradU(NdimE,NdimE) , DGradU(NdimE,NdimE), F(NdimE,NdimE)
		Real*8 :: D(NdimE,NdimE,NdimE,NdimE), Ppk(NdimE,NdimE)
		real*8 :: dVarPhi(NdimE), beta_Ten_dVarPhi(NdimE,NdimE), Ts(NdimE), Dts(NdimE,NdimE)
		real*8 :: Abu(NdimE,NodG*NdimE), Abb(NdimE,NdimE), Bb(NdimE), AbbInv(NdimE,NdimE), detAbb
		real*8 :: Bmat(NdimE*NdimE,NodG*NdimE), Dmat(NdimE*NdimE,NdimE*NdimE) , Pmat(NdimE*NdimE), & 
				  BmatT(NodG*NdimE,NdimE*NdimE), BmatBeta(NdimE*NdimE,NdimE), BmatBetaT(NdimE,NdimE*NdimE) 
		
		
!~ 		write(0,*) 'beta = ' , Beta 
!~ 		write(0,*) 'Param = ' , Param(1:7)
!~ 		write(0,*) 'matPar = ' , matPar
!~ 		write(0,*) 'damagePar = ' , damagePar
!~ 		write(0,*) 'SolU = ' , SolU
!~ 		write(0,*) 'SolDeltaU = ' , SolDeltaU
!~ 		write(0,*) 'cTs, aTs, damageTs, AreaS = ', cTs, aTs, damageTs, AreaS
!~ 		write(0,*) 'constLaw, iFtype,SolitaryNode, NodG, NGP, iDamageParType = ', constLaw, iFtype,SolitaryNode, NodG, NGP, iDamageParType
!~ 		pause	
		
		Ts = (1.0d0 - damageTs)*cTs*Beta
		Dts = (1.0d0 - damageTs)*cTs*deltaKron(1:NdimE,1:NdimE)
			
		do nG = 1, NGP ! LoopGauss
			
			call PtG%calcGradU(GradU,SolU - SolDeltaU,nG)
			
			call getF(F,iFtype) !!! with identity summed
			
			dVarPhi = PtG%dPhi_G(:,SolitaryNode,nG) !! solitary node is the last one, convention
			
			call VotimesW(beta_Ten_dVarPhi,beta,dVarPhi)  
			
			F = F + GradU - beta_Ten_dVarPhi
					
			call calcPpk(Ppk,F,NdimE,matPar,constLaw)		
			call calcD(D,F,NdimE,matPar,constLaw)	
			
			if(iDamageParType>0) then
				call damageUpdateEquilibrium(F,matPar,damagePar,Param,constLaw,nG) !~ 		
				call damageModifyTangent(D,Ppk,damagePar,Param, nG, 1)
			end if						
							
			call PtG%getBmat(Bmat,nG)
			BmatT = transpose(Bmat)				
			call voigtTen4toTen2(Dmat,D)
			call voigtTen2toVec(Pmat,Ppk)
			
			BmatBeta = Bmat(:,(solitaryNode-1)*NdimE + 1:(solitaryNode-1)*NdimE + 2) 
			BmatBetaT = transpose(BmatBeta) 
			
			dV=ptG%dV(nG)
									
			Abu = -matmul(BmatBetaT,matmul(Dmat,Bmat))*dV		
			Bb = matmul(BmatBetaT,Pmat)*dV - Ts*AreaS 
			Abb = matmul(BmatBetaT,matmul(Dmat,BmatBeta))*dV + Dts*AreaS
		
		end do !LoopGauss
		
		call matInv(AbbInv,detAbb,Abb)
		
		Beta = Beta + matmul(Abb,Bb - matmul(Abu,SolDeltaU))
!~ 		Beta = Beta + 0.0d0
		

	end subroutine

end module
