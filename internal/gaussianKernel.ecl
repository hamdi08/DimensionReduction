/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
/**
  *  Calculates the Gaussian kernel of the given matrix.
  *  @param X           Matrix X in DATASET(Layout_Cell) form.
  *  @param sigma       Sigma of the gaussian kernel in REAL.
  *  @return            Matrix in DATASET(Layout_Cell) form.
  *  @see               PBblas/Types.Layout_Cell
  */
  
 EXPORT gaussianKernel(DATASET(PBblas.Types.Layout_cell) X, REAL sigma=1.0) := FUNCTION
  X_dims := PBblas.internal.matdims.FromCells(X);
  nR := X_dims[1].m_rows;
  v1 := PBblas.HadamardProduct(X, X);
  v2 := TABLE(v1, {wi_id, x, sumR := SUM(GROUP, v)}, wi_id, x);
  v3 := PROJECT(v2, TRANSFORM(PBblas.Types.Layout_cell, 
                              SELF.y := 1,
                              SELF.v := LEFT.sumR,
                              SELF := LEFT));

 loopbody1(DATASET(PBblas.Types.Layout_Cell) m, UNSIGNED4 i) := FUNCTION
   m1 := PROJECT(m(y=i), TRANSFORM(PBblas.Types.Layout_cell,
                                  SELF.y:= i+1, SELF := LEFT));
   RETURN m+m1;
 END; 
 numIter := nR-1;
 v4 := LOOP(v3, COUNTER <= numIter, loopbody1(ROWS(LEFT), COUNTER));
 v4_5 := PBblas.tran(1.0, v3);
 loopbody2(DATASET(PBblas.Types.Layout_Cell) m, UNSIGNED4 i) := FUNCTION
   m1 := PROJECT(m(x=i), TRANSFORM(PBblas.Types.Layout_cell,
                                   SELF.x:= i+1, SELF := LEFT));
   RETURN m+m1;
 END; 
 v5 := LOOP(v4_5, COUNTER <= numIter, loopbody2(ROWS(LEFT), COUNTER));
 v6 := PBblas.axpy(1.0, v4, v5);
 v7 := PBblas.gemm(FALSE, TRUE, 1.0, X, X);
 v8 := PBblas.scal(-2, v7);
 v9 := PBblas.axpy(1.0, v6, v8);
 v10 := PROJECT(v9, TRANSFORM(PBblas.Types.Layout_cell,
                              SELF.v:= IF(LEFT.v < 0, 0, LEFT.v), SELF := LEFT));
 value_t := PBblas.Types.value_t;
 dimension_t := PBblas.Types.dimension_t;
 value_t expIt(value_t v, dimension_t r, dimension_t c) := EXP(v);
 v11 := PBblas.scal(-1 /(2 * sigma * sigma), v10);
 v12 := PBblas.Apply2Elements(v11, expIt);
 RETURN v12;
END;