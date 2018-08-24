/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
/**
  *  This module handles transition from Layout_cell format to Layout_cell_ext and vice versa.
  *  <p>Layout_cell is the matrix format of PBblas.
  *  <p>Layout_cell_ext is format that is used to store set of matrices.
  */
EXPORT PBMU := MODULE
 EXPORT To(DATASET(PBblas.Types.Layout_Cell) d, PBblas.Types.t_mu_no num) := 
  PROJECT(d, TRANSFORM(PBblas.Types.Layout_Cell_ext, SELF.mat_id := num, SELF := LEFT));
 EXPORT From(DATASET(PBblas.Types.Layout_Cell_ext) d, PBblas.Types.t_mu_no num) := 
  PROJECT(d(mat_id=num), TRANSFORM(PBblas.Types.Layout_Cell, SELF := LEFT));
END;