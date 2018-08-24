/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
IMPORT $ as dr;
/**
  *  Thin down a sparse matrix by removing any elements which are almost 0, while preserving the dimension 
  */
EXPORT thin(DATASET(PBblas.Types.Layout_cell) d) := FUNCTION
 dims := PBblas.internal.matdims.FromCells(d);
 thinD := d(ABS(v) > dr.config.RoundingError);
 RETURN dr.setDimension(thinD, dims[1].m_rows, dims[1].m_cols);
END;