/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
IMPORT $ AS dr;
/**
  *  Reduces the dimensionality of a matrix by projection on the linear/Gaussian/polynomial kernel-based principal components.
  *  @param X           Matrix X in DATASET(Layout_Cell) form.
  *  @param d           Number of dimensions to reduce to.
  *  @param tp          Type of the kernel in STRING: linear, Gaussian or polynomial; default 'linear'.
  *  @param param       Sigma for gaussian kernel or degree for polynomial kernel; default 1.
  *  @return            Matrix in DATASET(Layout_Cell) form.
  *  @see               PBblas/Types.Layout_Cell
  */
EXPORT kPCA(DATASET(PBblas.Types.Layout_Cell) X, UNSIGNED d=0, STRING tp='linear', REAL param=1) := FUNCTION
  value_t := PBblas.Types.value_t;
  dimension_t := PBblas.Types.dimension_t;
  X_dims := PBblas.internal.matdims.FromCells(X);
  m := X_dims[1].m_rows;
  n := X_dims[1].m_cols;
  
  K0 := dr.kernel(X, tp, param);
  
  //Centralizing K0
  onesM := dr.internal.ones(m,m);
  onesM_scl := PBblas.scal(1/m, onesM);
  t1 := PBblas.gemm(FALSE, FALSE, 1.0, K0, onesM_scl);
  t2 := PBblas.gemm(FALSE, FALSE, 1.0, onesM_scl, t1);
  t3 := PBblas.gemm(FALSE, FALSE, 1.0, onesM_scl, K0);
  t4 := PBblas.axpy(1.0, K0, T2);
  t5 := PBblas.axpy(1.0, T1, T3);
  K := PBblas.axpy(-1.0, T5, T4);
  K_scl := PBblas.scal(1/m, K);
  
  //Eigendecomposition of the scaled kernel matrix
  evecs := dr.eig(K_scl).vectors;
  evals_diag := dr.eig(K_scl).valuesM;
  evals := PBblas.vec.FromDiag(evals_diag);
  
  //sorting the eigenvectors according to the descending order of the eigenvalues
  Layout_Cell_eig_vec_val := RECORD(PBblas.Types.Layout_Cell)
   PBblas.Types.value_t eigval;
  END;
  evec_val := JOIN(evecs, evals, LEFT.wi_id = RIGHT.wi_id AND LEFT.y = RIGHT.x,
                   TRANSFORM(Layout_Cell_eig_vec_val,  
                             SELF.eigval := RIGHT.v, SELF := LEFT), LOOKUP);
  sorted_eig_vals_vecs := SORT(evec_val, -eigval, x);
  Layout_Cell_eig_vec_val_id := RECORD(Layout_Cell_eig_vec_val)
   UNSIGNED4 id; 
  END;
  Layout_cell_id := RECORD(PBblas.Types.Layout_cell)
   UNSIGNED id;
  END;
  evecs_id := PROJECT(evecs, TRANSFORM(Layout_cell_id, SELF.id := COUNTER, SELF := LEFT));
  sorted_eig_vals_vecs_id := PROJECT(sorted_eig_vals_vecs, 
                                     TRANSFORM(Layout_Cell_eig_vec_val_id, 
                                     SELF.id := COUNTER, SELF := LEFT));
  evecs_sorted := JOIN(evecs_id, sorted_eig_vals_vecs_id,
                       LEFT.wi_ID = RIGHT.wi_id AND LEFT.id = RIGHT.id,
                       TRANSFORM(PBblas.Types.Layout_cell,
                                 SELF.wi_id := LEFT.wi_id,
                                 SELF.x := LEFT.x,
                                 SELF.y := LEFT.y,  
                                 SELF.v := RIGHT.v), LOOKUP);
  
  //normalization of the eigenvectors
  value_t squareIt(value_t v, dimension_t r, dimension_t c) := v*v;
  evecs_sorted_sq := PBblas.Apply2Elements(evecs_sorted, squareIt);
  col_sum := TABLE(evecs_sorted_sq, {wi_id, y, sumC := SUM(GROUP, v)}, wi_id, y);
  evecs_sorted_sq_col_sum := PROJECT(col_sum, TRANSFORM(PBblas.Types.Layout_cell, 
                                                        SELF.x := 1,
                                                        SELF.v := LEFT.sumC,
                                                        SELF := LEFT));
  value_t sqrtIt(value_t v, dimension_t r, dimension_t c) := SQRT(v);
  norm_eigVector := PBblas.Apply2Elements(evecs_sorted_sq_col_sum, sqrtIt);
  dims_eigvec := PBblas.internal.matdims.FromCells(evecs_sorted);
  nR := dims_eigvec[1].m_rows;
  numIter := nR - 1;
  loopbody(DATASET(PBblas.Types.Layout_Cell) m, UNSIGNED4 i) := FUNCTION
    m1 := PROJECT(m(x=i), TRANSFORM(PBblas.Types.Layout_cell,
                                    SELF.x:= i+1, SELF := LEFT));
    RETURN m+m1;
  END; 
  rept_norm_eigVector := LOOP(norm_eigVector, COUNTER <= numIter, loopbody(ROWS(LEFT), COUNTER));
  value_t reciprocalIt(value_t v, dimension_t r, dimension_t c) := 1/v;
  rec_rept_norm_eigVector := PBblas.Apply2Elements(rept_norm_eigVector, reciprocalIt);
  norm_eigvector2 := PBblas.HadamardProduct(evecs_sorted, rec_rept_norm_eigVector);
  
  //getting only d dimensions
  kpc := norm_eigVector2(y<=d);
  Zspace := PBblas.gemm(FALSE, FALSE, 1.0, K0, kpc);
 
  RETURN Zspace;
END;