/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
IMPORT $ as dr;
/**
  *  Reduces the dimensionality of a matrix by projection on the principal components.
  *  @param A           Matrix A in DATASET(Layout_Cell) form.
  *  @param comp_cnt    Number of dimensions to reduce to.
  *  @return            Matrix in DATASET(Layout_Cell) form.
  *  @see               PBblas/Types.Layout_Cell
  */
EXPORT pca(DATASET(PBblas.Types.Layout_Cell) A, UNSIGNED comp_cnt=0) := MODULE
 CovA := dr.cov(A);
 U := dr.svd(CovA).UComp;
 // Ureduce - orthogonal vectors representing new space
 EXPORT Ureduce := IF(comp_cnt=0, U, U(y<=comp_cnt));
 //z-Normalizing matrix A
 Agg := TABLE(A, {wi_id, y, var := VARIANCE(GROUP, v), mean := AVE(GROUP, v)}, wi_id,y);
 A_norm := JOIN(A, Agg, 
                LEFT.wi_id=RIGHT.wi_id AND LEFT.y= RIGHT.y, 
                TRANSFORM(PBblas.Types.Layout_cell, 
                          SELF.v := (LEFT.v - RIGHT.mean)/SQRT(RIGHT.var), 
                          SELF := LEFT), 
                          LOOKUP);
 // ZComp - original features projected to the Ureduce space
 EXPORT ZComp := PBblas.gemm(FALSE, FALSE, 1.0, A_norm, Ureduce);
 //RETURN ZComp;
END;