
!> \ingroup h5block_data_f
!! See \ref H5Block3dWriteScalarFieldFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_write_scalar_field_r8 ( filehandle, name, buffer )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    REAL*8, INTENT(IN) :: buffer(*)    !< the array of data
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dReadScalarFieldFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_read_scalar_field_r8 ( filehandle, name, buffer )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    REAL*8, INTENT(OUT) :: buffer(*)   !< buffer to read the data into
END FUNCTION
 
!> \ingroup h5block_data_f
!! See \ref H5Block3dWriteVector3dFieldFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_write_vector3d_field_r8 ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    REAL*8, INTENT(IN) :: x(*)       !< the array of x data to write
    REAL*8, INTENT(IN) :: y(*)       !< the array of y data to write
    REAL*8, INTENT(IN) :: z(*)       !< the array of z data to write
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dReadVector3dFieldFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_read_vector3d_field_r8 ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    REAL*8, INTENT(OUT) :: x(*)      !< buffer to read the x data into
    REAL*8, INTENT(OUT) :: y(*)      !< buffer to read the y data into
    REAL*8, INTENT(OUT) :: z(*)      !< buffer to read the z data into
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dWriteScalarFieldFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_write_scalar_field_r4 ( filehandle, name, buffer )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    REAL*4, INTENT(IN) :: buffer(*)    !< the array of data
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dReadScalarFieldFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_read_scalar_field_r4 ( filehandle, name, buffer )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    REAL*4, INTENT(OUT) :: buffer(*)   !< buffer to read the data into
END FUNCTION
 
!> \ingroup h5block_data_f
!! See \ref H5Block3dWriteVector3dFieldFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_write_vector3d_field_r4 ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    REAL*4, INTENT(IN) :: x(*)       !< the array of x data to write
    REAL*4, INTENT(IN) :: y(*)       !< the array of y data to write
    REAL*4, INTENT(IN) :: z(*)       !< the array of z data to write
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dReadVector3dFieldFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_read_vector3d_field_r4 ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    REAL*4, INTENT(OUT) :: x(*)      !< buffer to read the x data into
    REAL*4, INTENT(OUT) :: y(*)      !< buffer to read the y data into
    REAL*4, INTENT(OUT) :: z(*)      !< buffer to read the z data into
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dWriteScalarFieldInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_write_scalar_field_i8 ( filehandle, name, buffer )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    INTEGER*8, INTENT(IN) :: buffer(*)    !< the array of data
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dReadScalarFieldInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_read_scalar_field_i8 ( filehandle, name, buffer )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    INTEGER*8, INTENT(OUT) :: buffer(*)   !< buffer to read the data into
END FUNCTION
 
!> \ingroup h5block_data_f
!! See \ref H5Block3dWriteVector3dFieldInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_write_vector3d_field_i8 ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    INTEGER*8, INTENT(IN) :: x(*)       !< the array of x data to write
    INTEGER*8, INTENT(IN) :: y(*)       !< the array of y data to write
    INTEGER*8, INTENT(IN) :: z(*)       !< the array of z data to write
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dReadVector3dFieldInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_read_vector3d_field_i8 ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    INTEGER*8, INTENT(OUT) :: x(*)      !< buffer to read the x data into
    INTEGER*8, INTENT(OUT) :: y(*)      !< buffer to read the y data into
    INTEGER*8, INTENT(OUT) :: z(*)      !< buffer to read the z data into
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dWriteScalarFieldInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_write_scalar_field_i4 ( filehandle, name, buffer )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    INTEGER*4, INTENT(IN) :: buffer(*)    !< the array of data
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dReadScalarFieldInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_read_scalar_field_i4 ( filehandle, name, buffer )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    INTEGER*4, INTENT(OUT) :: buffer(*)   !< buffer to read the data into
END FUNCTION
 
!> \ingroup h5block_data_f
!! See \ref H5Block3dWriteVector3dFieldInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_write_vector3d_field_i4 ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    INTEGER*4, INTENT(IN) :: x(*)       !< the array of x data to write
    INTEGER*4, INTENT(IN) :: y(*)       !< the array of y data to write
    INTEGER*4, INTENT(IN) :: z(*)       !< the array of z data to write
END FUNCTION

!> \ingroup h5block_data_f
!! See \ref H5Block3dReadVector3dFieldInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_3d_read_vector3d_field_i4 ( filehandle, name, x, y, z )
    INTEGER*8, INTENT(IN) :: filehandle  !< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: name !< the name of the dataset
    INTEGER*4, INTENT(OUT) :: x(*)      !< buffer to read the x data into
    INTEGER*4, INTENT(OUT) :: y(*)      !< buffer to read the y data into
    INTEGER*4, INTENT(OUT) :: z(*)      !< buffer to read the z data into
END FUNCTION

!> \ingroup h5block_attrib_f
!! See \ref H5BlockWriteFieldAttribFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_writefieldattrib_r8 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle		!< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name	!< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name	!< the name of the attribute
    REAL*8, INTENT(IN) :: buffer(*)		!< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelems		!< the number of elements in the array
END FUNCTION

!> \ingroup h5block_attrib_f
!! See \ref H5BlockReadFieldAttribFloat64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_readfieldattrib_r8 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle		!< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name	!< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name	!< the name of the attribute
    REAL*8, INTENT(IN) :: buffer(*)		!< the buffer to read into
END FUNCTION

!> \ingroup h5block_attrib_f
!! See \ref H5BlockWriteFieldAttribFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_writefieldattrib_r4 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle		!< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name	!< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name	!< the name of the attribute
    REAL*4, INTENT(IN) :: buffer(*)		!< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelems		!< the number of elements in the array
END FUNCTION

!> \ingroup h5block_attrib_f
!! See \ref H5BlockReadFieldAttribFloat32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_readfieldattrib_r4 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle		!< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name	!< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name	!< the name of the attribute
    REAL*4, INTENT(IN) :: buffer(*)		!< the buffer to read into
END FUNCTION

!> \ingroup h5block_attrib_f
!! See \ref H5BlockWriteFieldAttribInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_writefieldattrib_i8 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle		!< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name	!< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name	!< the name of the attribute
    INTEGER*8, INTENT(IN) :: buffer(*)		!< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelems		!< the number of elements in the array
END FUNCTION

!> \ingroup h5block_attrib_f
!! See \ref H5BlockReadFieldAttribInt64
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_readfieldattrib_i8 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle		!< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name	!< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name	!< the name of the attribute
    INTEGER*8, INTENT(IN) :: buffer(*)		!< the buffer to read into
END FUNCTION

!> \ingroup h5block_attrib_f
!! See \ref H5BlockWriteFieldAttribInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_writefieldattrib_i4 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle		!< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name	!< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name	!< the name of the attribute
    INTEGER*4, INTENT(IN) :: buffer(*)		!< the array of data to write into the attribute
    INTEGER*8, INTENT(IN) :: nelems		!< the number of elements in the array
END FUNCTION

!> \ingroup h5block_attrib_f
!! See \ref H5BlockReadFieldAttribInt32
!! \return 0 on success or error code
!<
INTEGER*8 FUNCTION h5bl_readfieldattrib_i4 ( filehandle, field_name, attrib_name, attrib_value, attrib_nelem)
    INTEGER*8, INTENT(IN) :: filehandle		!< the handle returned at file open
    CHARACTER(LEN=*), INTENT(IN) :: field_name	!< the name of the field
    CHARACTER(LEN=*), INTENT(IN) :: attrib_name	!< the name of the attribute
    INTEGER*4, INTENT(IN) :: buffer(*)		!< the buffer to read into
END FUNCTION
