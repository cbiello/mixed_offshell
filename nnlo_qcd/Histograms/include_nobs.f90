  integer, public, parameter :: Nobs = 22  ! total number of observabels MUST BE SMALLER THAN NOBS_MAX = 60
  integer, public, parameter :: Nnoncum = 21 ! total number of non cumulative histograms (must be in the first position)
  character :: cut_id*(50) = 'fid1'
