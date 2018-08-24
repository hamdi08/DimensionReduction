/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
IMPORT $ AS dr;
/**
  *  Decomposes a matrix of shape (m,n) by Singular Value Decomposition.
  *  <p>Implements SVD(A)=U*S*transpose(V).
  *  <p>A is a matrix with m rows and n columns.
  *
  *  @param A         Matrix A in DATASET(Layout_Cell) form with shape (m,n).
  *  @return UComp    Left singular matrix of shape (m,m).
  *  @return S2Comp   Diagonal matrix of eigenvalues of shape (m,n).
  *  @return SComp    Diagonal matrix of singular values of shape (m,n).
  *  @return VComp    Right singular matrix of shape (n,n).
  *  @see             PBblas/Types.Layout_Cell
  */
EXPORT svd(DATASET(PBblas.Types.Layout_cell) A) := MODULE
 SHARED AAt := PBblas.gemm(FALSE, TRUE, 1.0, A, A);
 EXPORT UComp := dr.eig(AAt).vectors;
 EXPORT S2Comp := dr.eig(AAt).valuesM;
 value_t := PBblas.Types.value_t;
 dimension_t := PBblas.Types.dimension_t;
 value_t squareRootIt(value_t v, dimension_t r, dimension_t c) := SQRT(v);
 EXPORT SComp := PBblas.Apply2Elements(S2Comp, squareRootIt);
 SHARED AtA := PBblas.gemm(TRUE, FALSE, 1.0, A, A);
 EXPORT VComp := dr.eig(AtA).vectors;
END;