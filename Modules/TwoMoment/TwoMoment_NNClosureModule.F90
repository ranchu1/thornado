MODULE TwoMoment_NNClosureModule

  USE KindModule, ONLY: DP, Pi, One, Two
  USE hdf5

  implicit none

  INTEGER, PARAMETER  :: input_n = 5
  INTEGER, PARAMETER  :: hid1 = 30
  INTEGER, PARAMETER  :: hid2 = 20
  INTEGER, PARAMETER  :: hid3 = 10
  INTEGER, PARAMETER  :: hid4 = 1

  INTEGER             :: hdferr
  INTEGER(HID_T)      :: file_id
  INTEGER(HID_T)      :: group_id1, group_id2

  real(dp) :: bias0(hid1)
  real(dp) :: bias1(hid2)
  real(dp) :: bias2(hid3)
  real(dp) :: bias3(hid4)

  real(dp) :: neuron0(hid1)
  real(dp) :: neuron1(hid2)
  real(dp) :: neuron2(hid3)
  real(dp) :: neuron3(hid4)

  real(dp) :: weights0(hid1,input_n)
  real(dp) :: weights1(hid2,hid1)
  real(dp) :: weights2(hid3,hid2)
  real(dp) :: weights3(hid4,hid3)

  integer(HSIZE_T), dimension(1) :: datesize_1d
  integer(HSIZE_T), dimension(2) :: datesize_2d

  contains

  subroutine nnclosure( inputs, k_NN )

    real(dp), intent(in)  :: inputs(input_n)
    real(dp), intent(out) :: k_NN

    neuron0 = forward( input_n, hid1, inputs,  weights0, bias0, 'tanh' )
    neuron1 = forward( hid1,    hid2, neuron0, weights1, bias1, 'tanh' )
    neuron2 = forward( hid2,    hid3, neuron1, weights2, bias2, 'tanh' )
    neuron3 = forward( hid3,    hid4, neuron2, weights3, bias3, 'sigmoid' )
    k_NN = neuron3(1)

  end subroutine


  subroutine InitializeNNClosure( FileName )

    character(*), intent(in) :: FileName

    write(*,*) 'InitializeNNClosure : ', FileName

    call OpenFileHDF( FileName, file_id )
  
    call OpenGroupHDF( 'dense', file_id, group_id1 )
    call OpenGroupHDF( 'dense', group_id1, group_id2 )
    datesize_1d = hid1
    call Read1dHDF_double( 'bias:0', bias0, group_id2, datesize_1d )
    datesize_2d(1) = hid1
    datesize_2d(2) = input_n
    call Read2dHDF_double( 'kernel:0', weights0, group_id2, datesize_2d )
    call CloseGroupHDF( group_id2 )
    call CloseGroupHDF( group_id1 )
  
    call OpenGroupHDF( 'dense_1', file_id, group_id1 )
    call OpenGroupHDF( 'dense_1', group_id1, group_id2 )
    datesize_1d = hid2
    call Read1dHDF_double( 'bias:0', bias1, group_id2, datesize_1d )
    datesize_2d(1) = hid2
    datesize_2d(2) = hid1
    call Read2dHDF_double( 'kernel:0', weights1, group_id2, datesize_2d )
    call CloseGroupHDF( group_id2 )
    call CloseGroupHDF( group_id1 )
  
    call OpenGroupHDF( 'dense_2', file_id, group_id1 )
    call OpenGroupHDF( 'dense_2', group_id1, group_id2 )
    datesize_1d = hid3
    call Read1dHDF_double( 'bias:0', bias2, group_id2, datesize_1d )
    datesize_2d(1) = hid3
    datesize_2d(2) = hid2 
    call Read2dHDF_double( 'kernel:0', weights2, group_id2, datesize_2d )
    call CloseGroupHDF( group_id2 )
    call CloseGroupHDF( group_id1 )

    call OpenGroupHDF( 'dense_3', file_id, group_id1 )
    call OpenGroupHDF( 'dense_3', group_id1, group_id2 )
    datesize_1d = hid3
    call Read1dHDF_double( 'bias:0', bias3, group_id2, datesize_1d )
    datesize_2d(1) = hid4
    datesize_2d(2) = hid3
    call Read2dHDF_double( 'kernel:0', weights3, group_id2, datesize_2d )
    call CloseGroupHDF( group_id2 )
    call CloseGroupHDF( group_id1 )
  
    call CloseFileHDF( file_id )

    write(*,*) 'hdferr:', hdferr ! 0 is correct

    if( hdferr == -1 )then
       write(*,*)
       write(*,*) 'Wrong in InitializeNNClosure: missing file '
       stop
    end if

  end subroutine InitializeNNClosure


  subroutine OpenFileHDF( FileName, file_id )

    character(len=*), intent(in) :: FileName
    INTEGER(HID_T), INTENT(out)  :: file_id

    call h5open_f(hdferr)
    call h5fopen_f( TRIM( FileName ), H5F_ACC_RDONLY_F, file_id, hdferr)

  end subroutine OpenFileHDF


  subroutine CloseFileHDF( file_id )

    INTEGER(HID_T), INTENT(in)   :: file_id

    call h5fclose_f( file_id, hdferr )

  end subroutine CloseFileHDF


  subroutine OpenGroupHDF( GroupName, file_id, group_id )

    CHARACTER(len=*), INTENT(in) :: GroupName
    INTEGER(HID_T), INTENT(in)   :: file_id
    INTEGER(HID_T), INTENT(out)  :: group_id

    call h5gopen_f( file_id, TRIM( GroupName ), group_id, hdferr )

  end subroutine OpenGroupHDF


  subroutine CloseGroupHDF( group_id )

    INTEGER(HID_T), INTENT(in) :: group_id

    call h5gclose_f( group_id, hdferr )

  end subroutine CloseGroupHDF


  subroutine Read1dHDF_double( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                    :: name
    INTEGER(HID_T)                               :: group_id
    INTEGER(HSIZE_T), DIMENSION(1), INTENT(in)   :: datasize
    REAL(dp), DIMENSION(:), INTENT(out)          :: values

    INTEGER(HID_T)                               :: dataset_id

    call h5dopen_f( group_id, name, dataset_id, hdferr )
    call h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, values, datasize, hdferr )
    call h5dclose_f( dataset_id, hdferr )

   end subroutine Read1dHDF_double


   subroutine Read2dHDF_double( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                    :: name
    INTEGER(HID_T)                               :: group_id
    INTEGER(HSIZE_T), DIMENSION(2), INTENT(in)   :: datasize
    REAL(dp), DIMENSION(:,:), INTENT(out)          :: values

    INTEGER(HID_T)                               :: dataset_id

    call h5dopen_f( group_id, name, dataset_id, hdferr )
    call h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, values, datasize, hdferr )
    call h5dclose_f( dataset_id, hdferr )

   end subroutine Read2dHDF_double


   function forward( Ni, No, inputs, weights, basis, activation )

     integer, intent(in)                    :: Ni, No ! number of in and out
     real(dp), dimension(Ni),    intent(in) :: inputs
     real(dp), dimension(No,Ni), intent(in) :: weights
     real(dp), dimension(Ni),    intent(in) :: basis
     character(*),               intent(in) :: activation
     
     real(dp), dimension(No)                :: forward

     integer :: i
     real(dp), dimension(No) :: Z
 
     Z = MATMUL(weights, inputs) + basis

     do i = 1, no
       if( trim( activation ) == 'sigmoid' ) then
         forward(i) = sigmoid( Z(i) )
       else if( trim( activation ) == 'tanh' ) then
         forward(i) = tanh( Z(i) )
       else
         stop 'unknown activation'
       end if
     end do

   end function forward


   function tanh( x )

     real(dp), intent(in) :: x
     real(dp)             :: tanh

     tanh = Two / ( One + fexp( - Two * x ) ) - One

   end function tanh


   function sigmoid( x )

     real(dp), intent(in) :: x
     real(dp)             :: sigmoid

     sigmoid = One / ( One + fexp( - x ) )

   end function sigmoid


   function fexp( x )

     real(dp), intent(in) :: x
     real(dp)             :: fexp

     real(dp), parameter  :: expmax = 300.d0
     real(dp), parameter  :: expmin = -300.d0

     fexp            = DEXP( DMIN1( expmax, DMAX1( expmin, x ) ) )

   end function fexp

END MODULE TwoMoment_NNClosureModule
