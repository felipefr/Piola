module fibresMod

use funcAux
use globalVariables, only : NdimE

implicit none

real*8 :: dummy
integer:: LengthParam,LengthJParam 
COMMON dummy,LengthParam,LengthJParam

integer , parameter :: 	Ipos_flag1 = 1 , Ipos_flag2 = 2	, Ipos_Lf = 3, &
							Ipos_Areaf = 4, Ipos_Vf = 5, Ipos_lfa = 6 , Ipos_af = 7, & 
							Ipos_yrel1 = 9, Ipos_yrel2 = 11
integer , parameter :: 	Ipos_measRVE = 1 , Ipos_measFibres = 2, Ipos_yG = 3 , Ipos_Bten = 5	, Ipos_BtenInvT = 9
integer :: i_fibresMod, j_fibresMod
real*8 :: yG(NdimE), measRVE=0.0d0, measFibres=0.0d0, Bten(NdimE,NdimE), Bten_invT(NdimE,NdimE)

DATA((Bten(i_fibresMod,j_fibresMod),i_fibresMod=1,NdimE),j_fibresMod=1,NdimE) /0.0d0,0.0d0,0.0d0,0.0d0/
DATA((Bten_invT(i_fibresMod,j_fibresMod),i_fibresMod=1,NdimE),j_fibresMod=1,NdimE) /0.0d0,0.0d0,0.0d0,0.0d0/
DATA(yG(i_fibresMod),i_fibresMod=1,NdimE) /0.0d0,0.0d0/

end module

