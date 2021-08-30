module interfaceFortran

	implicit none
	
	public executerElementC, executerSymbolicC
	private fortranMatrix2Carray, fortranMatrix2CarrayI
	
	contains

	subroutine executerElementC(id_Elem_Family, AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, & 
								Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)

		integer , intent(in) :: id_Elem_Family,MaxLRows,Ndim,iDofT, NodElt! all integers
		Real*8 , intent(in) :: DelT, DTm,Time ! all reals
		Real*8 , intent(out) :: AE(MaxLRows*MaxLRows), BE(MaxLRows)
		Real*8 , intent(in) :: XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)
		
		real*8 :: AEaux(MaxLRows,MaxLRows)
		
		AEaux = 0.0d0
		BE = 0.0d0
		call executerElement(id_Elem_Family, AEaux, BE, MaxLRows, XLL, NDim, iDofT, NodElt, & 
								Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
		
		
		call fortranMatrix2Carray(AE,AEaux,MaxLRows,MaxLRows)
		
	end subroutine

	subroutine executerSymbolicC(id_Elem_Family, Coupling,CommonPar,iDofT,Ndim,MaxLRows)

		Integer, intent(in) :: id_Elem_Family, MaxLRows,Ndim,iDofT
		Integer, intent(out) :: Coupling(MaxLRows,MaxLRows)
		Real*8 , intent(in) :: CommonPar(*)
		
		integer iAdd 
		integer :: CouplingAux(MaxLRows,MaxLRows)
		
		iAdd = 0;
		CouplingAux = 0;
		
		call executerSymbolic(id_Elem_Family, CouplingAux,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		
		call fortranMatrix2CarrayI(Coupling,CouplingAux,MaxLRows,MaxLRows)
		
	end subroutine 

	subroutine fortranMatrix2Carray(VinC,MinF,n,m)
		real*8 , intent(out) :: VinC(n*m)
		real*8 , intent(in) :: MinF(n,m)
		integer , intent(in) :: n,m
		integer i , j , ip
		
		do j  = 1,m
			do i  = 1,n
				ip = (i-1)*m + j
				VinC(ip) = MinF(i,j)
			end do
		end do
		
	end subroutine

	subroutine fortranMatrix2CarrayI(VinC,MinF,n,m)
		integer , intent(out) :: VinC(n*m)
		integer , intent(in) :: MinF(n,m)
		integer , intent(in) :: n,m
		integer i , j , ip
		
		do j  = 1,m
			do i  = 1,n
				ip = (i-1)*m + j
				VinC(ip) = MinF(i,j)
			end do
		end do
		
	end subroutine
	
end module

