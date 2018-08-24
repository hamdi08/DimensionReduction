/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
IMPORT $ as dr;
/**
  *  Decomposes a  real square matrix of shape (n,n) by QR decomposition.
  *  <p>Implements QR(A)=Q*R.
  *  <p>A is a square matrix with n rows and n columns.
  *  <p>Q is a orthogonal matrix of shape (n,n).
  *  <p>R is a upper triangular matrix of shape (n,n).
  *
  *  @param matrix          Matrix in DATASET(Layout_Cell) form with shape (n,n).
  *  @return [R,Q]          Set of two matrices R and Q.
  *  @see                   PBblas/Types.Layout_Cell
  */
qr_comp := ENUM(Q = 1, R = 2);
EXPORT DATASET(PBblas.Types.Layout_Cell_ext) QR(DATASET(PBblas.Types.Layout_Cell) matrix) := FUNCTION
 dims := PBblas.internal.matdims.FromCells(matrix);
 n := dims[1].m_rows;
 loopBody(DATASET(PBblas.Types.Layout_Cell_ext) ds, UNSIGNED4 k) := FUNCTION
  Q := dr.internal.PBMU.From(ds, qr_comp.Q);
  R := dr.internal.PBMU.from(ds, qr_comp.R);
  hModule := dr.internal.Householder(PBblas.vec.FromCol(R, k), k, n);
  hM := hModule.Reflection();
  R1 := dr.internal.PBMU.To(PBblas.gemm(FALSE, FALSE, 1.0, hM, R), qr_comp.R);
  Q1 := dr.internal.PBMU.To(PBblas.gemm(FALSE, FALSE, 1.0, Q, hM), qr_comp.Q);
  RETURN R1+Q1;
 END;
 RETURN LOOP(dr.internal.PBMU.To(matrix, qr_comp.R) + 
             dr.internal.PBMU.To(PBblas.eye(n), qr_comp.Q), 
             n-1, 
             loopBody(ROWS(LEFT),COUNTER));
END;