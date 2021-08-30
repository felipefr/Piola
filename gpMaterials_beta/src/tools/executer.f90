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
		Case (502)
			Call enforcePeriodic2D( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (506)
			Call computeTotalDisp2D( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (508)
			Call PostProcessingElem2D( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (511)
			Call computeTangHom( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (512)
			Call minRestrictionBC2DExact( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
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
		Case (520)
			Call DamageGeneric( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (522)
			Call finiteStrainNew( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (523)
			Call flucTangHom( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (524)
			Call postProcessingRVEboundary( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
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
		Case (540)
			Call computePressure( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (541)
			Call flucTangHomGen( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (542)
			Call computeTangHomGen( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (600)
			Call DecisionBifurcation( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (601)
			Call LocDomainBC( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (700)
			Call computeMinDetQ( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (701)
			Call MarkLocPoint( AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (702)
			Call LocPointsBC(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (797)
			Call posProcInterp(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
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
		Case (802)
			Call TotalDisp(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (803)
			Call DamageGenericGen(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (804)
			Call canonicalProblem(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (805)
			Call TangentHom(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (806)
			Call PosProcElem(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (807)
			Call minRestrictionRVE(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (900)
			Call insertDeformation(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (901)
			Call setMaterialParam(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (902)
			Call setDamageParam(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (903)
			Call insertDeformation2(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (904)
			Call insertDeformationSimple(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time) 
		Case (1001)
			Call hyperFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1002)
			Call fibresTraction(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1003)
			Call fibresHom(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1004)
			Call fibresDisp(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1005)
			Call fibresHom2(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1006)
			Call linearFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1007)
			Call reactionLag(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1008)
			Call NodalForce(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1009)
			Call BarElements(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1010)
			Call nonlinearFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1011)
			Call linearFibres2(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1012)
			Call nonlinearFibres_qfBased(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1013)
			Call damageGenericFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1014)
			Call minRestrictionFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1015)
			Call computeAnisotropyTensor(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1016)
			Call computeAnisotropyTensorInv(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1017)
			Call minRestrictionFibresAniso(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1018)
			Call computeAnisoTensors(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1019)
			Call computeAnisoTensorsInv(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1020)
			Call RestrictionsFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1021)
			Call fibresHom3(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1022)
			Call posProcFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1023)
			Call MRplusZeroMean(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1024)
			Call StressHomGen(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1025)
			Call zeroAverageFibre(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1026)
			Call linearBarElementMS(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1027)
			Call fibresHom4(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1028)
			Call crossLinks_restriction(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)		
		Case (1029)
			Call computeAnisotropyTensor_boundary(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (1030)
			Call minRestrictionFibresAniso_boundary(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1031)
			Call minRestrictionContinuumFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1032)
			Call initParamFibres(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1033)
			Call nonlinearFibresVoigt(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (1034)
			Call networkConstraint(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)

		Case (2000)
			Call AdjointProblem(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (2001)
			Call identifyUpdateDamage(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (2002)
			Call computeSens(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (2003)
			Call fillParamSensibility(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)						
		Case (2004)
			Call fillParamSensibility_GradD(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (2005)
			Call mechResidual(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (2006)
			Call fillParamSensibility_mechRes(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)											
		Case (3000)
			Call FSgen_EFEM(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (3001)
			Call FSgen_EFEM_voigt(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (3002)
			Call FSgen_EFEM_condesated(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (3003)
			Call LinElas_EFEM_voigt(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (3004)
			Call LinElas_EFEM_condesated(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)	
		Case (3005)
			Call LinElas_beta_update(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)
		Case (3006)
			Call damageGenericGenLinElas(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
								Sol0,Sol1,CommonPar,&
								Param,JParam,DelT,DTm,Time)									
		Case (3007)
			Call computeMinDetQLinElas(AE,BE,MaxLRows,XLL,NDim,iDofT,NodElt,&
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
		Case (502)
		!> 502 - enforcePeriodic2D.f90
			Call enforcePeriodic2DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (506)
		!> 506 - computeTotalDisp2D.f90
			Call computeTotalDisp2DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (508)
		!> 508 - PostProcessingElem2D.f90
			Call PostProcessingElem2DS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (511)
		!> 511 - computeTangHom.f90
			Call computeTangHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (512)
		!> 512 - minRestrictionBC2DExact.f90
			Call minRestrictionBC2DExactS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (518)
		!> 518 - DirichletNode.f90
			Call DirichletNodeS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (519)
		!> 519 - DirichletNodeLag.f90
			Call DirichletNodeLagS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (520)
		!> 520 - DamageGeneric.f90
			Call DamageGenericS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (522)
		!> 522 - finiteStrainNew.f90
			Call finiteStrainNewS (Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (523)
		!> 523 - flucTangHom.f90
			Call flucTangHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (524)
			Call postProcessingRVEboundaryS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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
		Case (540)
			Call computePressureS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (541)
			Call flucTangHomGenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (542)
			Call computeTangHomGenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (600)
			Call DecisionBifurcationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (601)
			Call LocDomainBCS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (700)
			Call computeMinDetQS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (701)
			Call MarkLocPointS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (702)
			Call LocPointsBCS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (797)
		!> 797 - posProcInterp.f90
			Call posProcInterpS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
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
		Case (802)
		!> 802 - TotalDisp.f90
			Call TotalDispS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! with inserted deformation
		Case (803)
		!> 803 - DamageGenericGen.f90
			Call DamageGenericGenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! with Material global
		Case (804)
		!> 804 - canonicalProblem.f90
			Call canonicalProblemS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! with Material global
		Case (805)
		!> 805 - TangentHom.f90
			Call TangentHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! with Material global
		Case (806)
		!> 806 - posProcElem.f90
			Call PosProcElemS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! with Material global
		Case (807)
		!> 807 - minRestrictionRVE.f90
			Call minRestrictionRVES(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! with Material global
		Case (900)
		!> 900 - insertDeformation.f90
			Call insertDeformationS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! set F properly
		Case (901)
		!> 901 - setMaterialParam.f90
			Call setMaterialParamS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! set material constitutive parameters once
		Case (902)
		!> 902 - setDamageParam.f90
			Call setDamageParamS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! set material constitutive parameters once
		Case (903)
		!> 903 - insertDeformation2.f90
			Call insertDeformation2S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! set material constitutive parameters once
		Case (904)
		!> 904 - insertDeformationSimple.f90
			Call insertDeformationSimpleS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd) !!! set material constitutive parameters once
		Case (1001)
		!> 1001 -  hyperFibres.f90
			Call hyperFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1002)
		!> 1002 - fibresTraction.f90
			Call fibresTractionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1003)
		!> 1003 - fibresHom.f90
			Call fibresHomS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1004)
		!> 1004 - fibresDisp.f90
			Call fibresDispS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1005)
		!> 1005 - fibresHom2.f90
			Call fibresHom2S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1006)
		!> 1006 - linearFibres.f90
			Call linearFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1007)
		!> 1007 - reactionLag.f90
			Call reactionLagS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1008)
		!> 1008 - NodalForce.f90
			Call NodalForceS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1009)
			Call BarElementsS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1010)
			Call nonlinearFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1011)
			Call linearFibres2S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1012)
		!> 1012 - nonlinearFibres_qfBased.f90
			Call nonlinearFibres_qfBasedS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1013)
		!> 1013 - damageGenericFibres.f90
			Call damageGenericFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1014)
		!> 1014 - minRestrictionFibres.f90
			Call minRestrictionFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1015)
		!> 1015 - computeAnisotropyTensor.f90
			Call computeAnisotropyTensorS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1016)
		!> 1016 - computeAnisotropyTensorInv.f90
			Call computeAnisotropyTensorInvS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1017)
		!> 1017 - minRestrictionFibresAniso.f90
			Call minRestrictionFibresAnisoS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1018)
		!> 1018 - computeAnisoTensors.f90
			Call computeAnisoTensorsS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1019)
		!> 1019 - computeAnisoTensorsInv.f90
			Call computeAnisoTensorsInvS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1020)
		!> 1020 - RestrictionsFibres.f90
			Call RestrictionsFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1021)
		!> 1021 - fibresHom3.f90
			Call fibresHom3S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)		
		Case (1022)
		!> 1022 - PosProcFibres.f90
			Call PosProcFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1023)
		!> 1023 - MRplusZeroMean.f90
			Call MRplusZeroMeanS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1024)
		!> 1024 - StressHomGen.f90
			Call StressHomGenS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1025)
		!> 1025 - zeroAverageFibre.f90
			Call zeroAverageFibreS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)		
		Case (1026)
		!> 1026 - linearBarElementMS.f90
			Call linearBarElementMSS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1027)
		!> 1027 - fibresHom4.f90
			Call fibresHom4S(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1028)
		!> 1028 - crossLinks_restriction.f90
			Call crossLinks_restrictionS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1029)
		!> 1029 - computeAnisotropyTensor_boundary.f90
			Call computeAnisotropyTensor_boundaryS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1030)
		!> 1030 - minRestrictionFibresAniso_Boundary.f90
			Call minRestrictionFibresAniso_boundaryS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1031)
		!> 1031 - minRestrictionContinuumFibres.f90
			Call minRestrictionContinuumFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1032)
		!> 1032 - initParamFibres.f90
			Call initParamFibresS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)	
		Case (1033)
		!> 1033 - nonlinearFibresVoigt.f90
			Call nonlinearFibresVoigtS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (1034)
		!> 1034 - networkConstraint.f90
			Call networkConstraintS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (2000)
		!> 2000 - AdjointProblem.f90
			Call AdjointProblemS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (2001)
		!> 2001 - identifyUpdateDamage.f90
			Call identifyUpdateDamageS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (2002)
		!> 2002 - computeSens.f90
			Call computeSensS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (2003)
		!> 2003 - fillParamSensibility.f90
			Call fillParamSensibilityS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (2004)
		!> 2004 - fillParamSensibility_GradD.f90
			Call fillParamSensibility_GradDS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (2005)
		!> 2005 - mechResidual.f90
			Call mechResidualS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (2006)
		!> 2006 - fillParamSensibility_mechRes.f90
			Call fillParamSensibility_mechResS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (3000)
		!> 3000 - FSgen_EFEM.f90
			Call FSgen_EFEMS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (3001)
		!> 3001 - FSgen_EFEM_voigt.f90
			Call FSgen_EFEM_voigtS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (3002)
		!> 3002 - FSgen_EFEM_condesated.f90
			Call FSgen_EFEM_condesatedS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (3003)
		!> 3003 - LinElas_EFEM_voigt.f90
			Call LinElas_EFEM_voigtS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (3004)
		!> 3004 - LinElas_EFEM_condesated.f90
			Call LinElas_EFEM_condesatedS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (3005)
		!> 3005 - LinElas_beta_update.f90
			Call LinElas_EFEM_condesatedS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (3006)
		!> 3006 - DamageGenericGenLinElas.f90
			Call DamageGenericGenLinElasS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
		Case (3007)
		!> 3007 - ComputeMinDetQLinElas.f90
			Call computeMinDetQLinElasS(Coupling,CommonPar,iDofT,Ndim,MaxLRows,iAdd)
    Case Default
           Write(*,*) '!!!!!!! WARNING, ELEMENT NOT IDENTIFIED'	, Id_Elem_Family
           Write(*,*) 'while trying to determine the Coupling Matrix'	


	End Select

end subroutine  



