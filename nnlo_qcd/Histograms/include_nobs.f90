  integer, public, parameter :: Nobs = 17  ! total number of observabels MUST BE SMALLER THAN NOBS_MAX = 60
  integer, public, parameter :: Nnoncum = 16 ! total number of non cumulative histograms (must be in the first position)
  character :: cut_id*(50) = 'fid1'
