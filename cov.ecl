/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
/**
  *  Calculates the covariance matrix of a given matrix after z-Normalization.
  *  <p>Implements (1/m)*transpose(A_norm)*A_norm.
  *  <p>A is a matrix with m rows and n columns.
  *  <p>A_norm is the z-normalized A with same shape.
  *
  *  @param A             Matrix A in DATASET(Layout_Cell) form with shape (m,n).
  *  @return Cov(A_norm)  Matrix in DATASET(Layout_Cell) form with shape (n,n).
  *  @see                 PBblas/Types.Layout_Cell
  */
EXPORT cov(DATASET(PBblas.Types.Layout_cell) A) :=  FUNCTION
 //z-Normalizing
 Agg := TABLE(A, {wi_id, y, var := VARIANCE(GROUP, v), mean := AVE(GROUP, v)}, wi_id, y);
 A_norm := JOIN(A, Agg, 
                LEFT.wi_id=RIGHT.wi_id AND LEFT.y= RIGHT.y, 
                TRANSFORM(PBblas.Types.Layout_cell, 
                          SELF.v := (LEFT.v - RIGHT.mean)/SQRT(RIGHT.var), 
                          SELF := LEFT), 
                          LOOKUP);
 A_dims := PBblas.internal.matdims.FromCells(A_norm);
 m := A_dims[1].m_rows;
 //n := X_dims[1].m_cols;
 //Covariance matrix calculation
 covmat_nscl := PBblas.gemm(TRUE, FALSE, 1.0, A_norm, A_norm);
 C := PBblas.scal(1/m, covmat_nscl);
 RETURN C;
END;