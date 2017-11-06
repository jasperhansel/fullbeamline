	SUBROUTINE four2(data,isign)
	USE nrtype
	USE nr, ONLY : fourrow
	IMPLICIT NONE
	COMPLEX(dpc), DIMENSION(:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(dpc), DIMENSION(size(data,2),size(data,1)) :: temp
	call fourrow(data,isign)
	temp=transpose(data)
	call fourrow(temp,isign)
	data=transpose(temp)
	END SUBROUTINE four2
