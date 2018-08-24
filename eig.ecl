/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
IMPORT $ as dr;
/**
  *  Decomposes a  real symmetric matrix of shape (n,n) by Eigendecomposition.
  *  <p>Implements Eig(A)=U*S*transpose(U).
  *  <p>A is a  square matrix with n rows and n columns.
  *
  *  @param A               Matrix A in DATASET(Layout_Cell) form with shape (n,n).
  *  @return valuesM        Diagonal matrix of Eigenvalues of shape (n,n).
  *  @return vectors        Matrix where the columns are the Eigenvectors (n,n).
  *  @return convergence    Value after convergence loop.
  *  @see                   PBblas/Types.Layout_Cell
  */

EXPORT eig(DATASET(PBblas.Types.Layout_cell) A, UNSIGNED4 iter=200) := MODULE
 SHARED eig_comp := ENUM (T = 1, Q = 2, conv = 3);
 SHARED qr_comp_id := ENUM (Q = 1,  R = 2);
 EXPORT DATASET(PBblas.Types.Layout_Cell_ext) QRalgorithm() := FUNCTION
  QR0 := dr.QR(A);
  Q0 := dr.internal.PBMU.From(QR0, qr_comp_id.Q);
  R0 := dr.internal.PBMU.From(QR0, qr_comp_id.R);
  T0 := PBblas.gemm(FALSE, FALSE, 1.0, R0, Q0);
  Conv0 := DATASET([{1,1,1,0}],PBblas.Types.Layout_cell);
  loopBody(DATASET(PBblas.Types.Layout_Cell_ext) ds, UNSIGNED4 k) := FUNCTION
   T := dr.internal.PBMU.From(ds, eig_comp.T);
   Q := dr.internal.PBMU.From(ds, eig_comp.Q);
   Conv := dr.internal.PBMU.From(ds, eig_comp.conv);
   bConverged := PBblas.vec.Norm(PBblas.vec.FromDiag(T,-1)) < dr.internal.config.RoundingError;
   QR1 := dr.QR(T);
   QComp := dr.internal.PBMU.From(QR1, qr_comp_id.Q);
   Q1 := dr.internal.thin(PBblas.gemm(FALSE, FALSE, 1.0, Q, QComp));
   RComp := dr.internal.PBMU.From(QR1, qr_comp_id.R);
   T1 := dr.internal.thin(PBblas.gemm(FALSE, FALSE, 1.0, RComp, QComp));
   Conv1 :=  PROJECT(Conv, TRANSFORM(PBblas.Types.Layout_cell, SELF.v:=k, SELF:=LEFT));
   RETURN IF(bConverged, 
              ds, 
              dr.internal.PBMU.To(T1, eig_comp.T) + 
              dr.internal.PBMU.To(Q1, eig_comp.Q)+ 
              dr.internal.PBMU.To(Conv1, eig_comp.conv));
  END;
  RETURN LOOP(dr.internal.PBMU.To(T0, eig_comp.T) + 
              dr.internal.PBMU.To(Q0, eig_comp.Q) + 
              dr.internal.PBMU.To(Conv0, eig_comp.conv), 
              iter, 
              loopBody(ROWS(LEFT),COUNTER));
 END;
EXPORT valuesM := dr.internal.PBMU.From(QRalgorithm(), eig_comp.T);
//EXPORT valuesV := PBblas.vec.FromDiag(valuesM());
EXPORT vectors := dr.internal.PBMU.From(QRalgorithm(), eig_comp.Q);
EXPORT convergence := dr.internal.PBMU.From(QRalgorithm(), eig_comp.conv)[1].v;
END;