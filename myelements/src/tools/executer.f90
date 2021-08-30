subroutine executerElement(id_Elem_Family, AE, BE, MaxLRows, XLL, NDim, iDofT, NodElt, & 
							Sol0, Sol1, CommonPar, Param, JParam, DelT, DTm, Time)
	implicit none
	
	integer :: id_Elem_Family,MaxLRows,Ndim,iDofT, NodElt! all integers
	Real*8 :: DelT, DTm,Time ! all reals
	Real*8 :: AE(MaxLRows,MaxLRows), BE(MaxLRows), XLL(*), Sol0(*), Sol1(*), Param(*), JParam(*), CommonPar(*)
	

	
	Select Case (Id_Elem_Family)
		Case (100)
			Call trivialElement( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (101)
			Call PolyElement( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (102)
			Call NullElement( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (103)
			Call Cell2Point( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (104)
			Call IncrementVariables( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (518)
			Call DirichletNode( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (519)
			Call DirichletNodeLag( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (522)
			Call finiteStrainNew( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (525)
			Call finiteStrainDelta( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (529)
			Call finiteStrainGen( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (530)
			Call FbarReduced( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (531)
			Call FiniteStrainSpatial( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (532)
			Call FbarMaterial( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (533)
			Call NeumannFiniteStrain( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (534)
			Call NeumannFiniteStrainSpatial( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (535)
			Call NeumannRefTraction( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (538)
			Call FbarSpatial( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (539)
			Call NodalForceInLine( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (798)
			Call FSgenMixed_voigt(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (799)
			Call FSgen_voigt(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (800)
			Call FSgen(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (801)
			Call FSFbar(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (901)
			Call setMaterialParam(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)  																			
		Case Default
          Write(*,*) '!!!!!!! WARNING, ELEMENT NOT IDENTIFIED', Id_Elem_Family

       End Select 

end subroutine

subroutine executerSymbolic(id_Elem_Family, Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
	implicit none
	
	Integer :: id_Elem_Family, MaxLRows,Ndim,iDofT, iAdd
    Integer   Coupling(MaxLRows,MaxLRows)
    Real*8 CommonPar(*)
	
	Select Case (Id_Elem_Family)

		Case (100)
		!> 100 - TrivialElement.f90
			Call  trivialElementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (101)
		!> 101 - PolyElement.f90
			Call  PolyElementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (102)
		!> 102 - NullElement.f90
			Call  NullElementS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (103)
		!> 103 - Cell2Point.f90
			Call  Cell2PointS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (104)
		!> 104 - IncrementVariables.f90
			Call  IncrementVariablesS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (518)
		!> 518 - DirichletNode.f90
			Call DirichletNodeS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (519)
		!> 519 - DirichletNodeLag.f90
			Call DirichletNodeLagS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (522)
		!> 522 - finiteStrainNew.f90
			Call finiteStrainNewS (Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (525)
			Call finiteStrainDeltaS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (529)
			Call finiteStrainGenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (530)
			Call FbarReducedS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (531)
			Call FiniteStrainSpatialS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (532)
			Call FbarMaterialS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (533)
		!> 533 - NeumannFiniteStrain.f90
			Call NeumannFiniteStrainS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (534)
		!> 534 - NeumannFiniteStrainSpatial.f90
			Call NeumannFiniteStrainSpatialS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (535)
		!> 535 - NeumannRefTraction.f90
			Call NeumannRefTractionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (538)
			Call FbarSpatialS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (539)
		!> 539 - NodalForceInLine.f90 !!! For Cook Problem load
			Call NodalForceInLineS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (798)
		!> 798 - FSgenMixed_voigt.f90
			Call FSgenMixed_voigtS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (799)
		!> 799 - FSgen_voigt.f90
			Call FSgen_voigtS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (800)
		!> 800 - FSgen.f90
			Call FSgenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! with inserted deformation
		Case (801)
		!> 801 - FSFbar.f90
			Call FSFbarS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! with inserted deformation
		Case (901)
		!> 901 - setMaterialParam.f90
			Call setMaterialParamS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! set material constitutive parameters once

		
    Case Default
           Write(*,*) '!!!!!!! WARNING, ELEMENT NOT IDENTIFIED'	, Id_Elem_Family
           Write(*,*) 'while trying to determine the Coupling Matrix'	


	End Select

end subroutine  



