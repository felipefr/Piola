subroutine homogenisationFibres(PK,ParamFile, iShift_volfib, iShift_afib, iShift_sfib, &
								Nelem, NParamCell, NdimE, NTimeSteps)
		
	use converterFuncs
	
	real*8 , intent(inout) :: PK(NdimE, NdimE, NTimeSteps)
	integer, intent(in) :: Nelem, NTimeSteps, NParamCell, NdimE, iShift_volfib , iShift_afib, iShift_sfib 
	character(len=40), intent(in) :: ParamFile
	
	call homogenisationFibres_(PK, ParamFile, Nelem, NTimeSteps, NParamCell, NdimE, iShift_volfib , iShift_afib, iShift_sfib)

end subroutine
