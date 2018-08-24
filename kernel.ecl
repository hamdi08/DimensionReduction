/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
IMPORT $ AS dr;
/**
  *  Directs to the corresponding kernel function according to the given argument.
  *  @param X           Matrix X in DATASET(Layout_Cell) form.
  *  @param d           Number of dimensions to reduce to.
  *  @param tp          Type of the kernel in STRING: linear, gaussian or polynomial; default 'linear'.
  *  @param param       Sigma for Gaussian kernel or degree for polynomial kernel; default 1.
  *  @return            Matrix in DATASET(Layout_Cell) form.
  *  @see               PBblas/Types.Layout_Cell
  */
EXPORT kernel(DATASET(PBblas.Types.Layout_Cell) X, STRING tp='linear', REAL param=1) := FUNCTION
 K := CASE(tp, 'linear' => PBblas.gemm(FALSE, TRUE, 1.0, X, X),
               'poly' => dr.internal.polyKernel(X, param),
               'gaussian' => dr.internal.gaussianKernel(X, param));
RETURN K;
END;