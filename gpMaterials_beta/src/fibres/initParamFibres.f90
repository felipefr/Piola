Subroutine initParamFibres(AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
!     ------------------------------------------------------------------
	use funcAux
	use globalVariables, only : NdimE
	use fibresMod
	
	implicit none

	!   ===== SUBROUTINE ARGUMENTS  =======
	integer :: MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)

!~ 	!   =====   END ARGUMENTS  =======
         
	Real*8 , parameter :: tol = 1.0e-8
	Real*8 :: af(NdimE), Lf, Vf, Areaf, lfa, yrel1(NdimE), yrel2(NdimE), yrel_otimes_af(NdimE,NdimE)
	real*8 , allocatable ::   Xel(:) !! all have dimension NodG*NdimE
	integer :: NodG
	integer :: flag1, flag2, op 

	op = int(CommonPar(1))
	
	NodG = 2	
	
	call getSliceAllocate(Xel,XLL,1,NodG,1 ,NdimE, Ndim)
		
		
	if(op == 1) then !! initialize first variables
		
		Areaf = CommonPar(2)
		lfa = CommonPar(3)
	
		af(1) = Xel(3) - Xel(1)
		af(2) = Xel(4) - Xel(2)  
		
		Lf = norm2(af)
		af = af/Lf
		
		Vf = Areaf*Lf
		
		flag1 = 0
		flag2 = 0
		
		if( dabs(xel(1) - 0.0d0)  < tol .or. dabs(xel(1) - 1.0d0)  < tol) then
			flag1 = 1
		else if ( dabs(xel(2) - 0.0d0)  < tol .or. dabs(xel(2) - 1.0d0)  < tol) then
			flag1 = 1
		end if

		if( dabs(xel(3) - 0.0d0)  < tol .or. dabs(xel(3) - 1.0d0)  < tol) then
			flag2 = 1
		else if ( dabs(xel(4) - 0.0d0)  < tol .or. dabs(xel(4) - 1.0d0)  < tol) then
			flag2 = 1
		end if
		
		measFibres = measFibres + Vf
		yG = yG + 0.5d0*Vf*( Xel(1:2) + Xel(3:4) )
		
		Param(LengthParam  + Ipos_flag1 ) = flag1
		Param(LengthParam  + Ipos_flag2 ) = flag2	
		Param(LengthParam  + Ipos_Lf ) = Lf
		Param(LengthParam  + Ipos_Areaf ) = Areaf
		Param(LengthParam  + Ipos_Vf ) = Vf
		Param(LengthParam  + Ipos_lfa ) = lfa
		Param(LengthParam  + Ipos_af ) = af(1)
		Param(LengthParam  + Ipos_af + 1 ) = af(2)
		
	else if(op==2) then !! compute yG (take average)
		yG = yG/measFibres
		
	else if(op==3) then !! compute yrel and B
				
 		flag1 = nint(Param(LengthParam + Ipos_flag1))
 		flag2 = nint(Param(LengthParam + Ipos_flag2))
		Vf = Param(LengthParam +  Ipos_Vf)
		Areaf = Param(LengthParam + Ipos_Areaf)
		af = Param(LengthParam + Ipos_af: LengthJParam + Ipos_af+1)
		
		yrel1 = Xel(1:2) - yG
		yrel2 = Xel(3:4) - yG
		
 	    if(flag1 == 1) then 		
			call VotimesW(yrel_otimes_af,yrel1,af)
			Bten = Bten - (Areaf/measFibres)*yrel_otimes_af
		end if

	 	if(flag2 == 1) then 		
			call VotimesW(yrel_otimes_af,yrel2,af)
			Bten = Bten + (Areaf/measFibres)*yrel_otimes_af
		end if
 		
		Param(LengthParam + Ipos_yrel1) = yrel1(1)
		Param(LengthParam + Ipos_yrel1 + 1) = yrel1(2)
		Param(LengthParam + Ipos_yrel2) = yrel2(1)
		Param(LengthParam + Ipos_yrel2 + 1) = yrel2(2)
	
	else if(op==4) then !! compute Bten_invT and set then to global param
		measRVE = commonPar(2)
		
		call MatInv(Bten_invT,dummy,Bten)
		
		Bten_invT = transpose(Bten_invT)
		
!~ 		Param(LengthParam + Ipos_measRVE) = measRVE
!~ 		Param(LengthParam + Ipos_measFibres) = measFibres
!~ 		Param(LengthParam + Ipos_yG) = yG(1)
!~ 		Param(LengthParam + Ipos_yG + 1) = yG(2)
!~ 		Param(LengthParam + Ipos_Bten) = Bten(1,1)
!~ 		Param(LengthParam + Ipos_Bten + 1) = Bten(1,2)
!~ 		Param(LengthParam + Ipos_Bten + 2) = Bten(2,1)
!~ 		Param(LengthParam + Ipos_Bten + 3) = Bten(2,2)
!~ 		Param(LengthParam + Ipos_BtenInvT) = Bten_invT(1,1)
!~ 		Param(LengthParam + Ipos_BtenInvT + 1) = Bten_invT(1,2)
!~ 		Param(LengthParam + Ipos_BtenInvT + 2) = Bten_invT(2,1)
!~ 		Param(LengthParam + Ipos_BtenInvT + 3) = Bten_invT(2,2)
	
	end if
	
	AE(1,1) = 1.0d0
	BE(1) = Sol1(1)
	
end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Subroutine initParamFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use funcAux
    
    IMPLICIT NONE
    Integer :: MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
    
    Coupling(1,1) = 1

end Subroutine
