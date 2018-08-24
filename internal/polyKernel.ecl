/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
/**
  *  Calculates the polynomial kernel of the given matrix.
  *  @param X           Matrix X in DATASET(Layout_Cell) form.
  *  @param deg         Degree of the polynomial kernel in REAL.
  *  @return            Matrix in DATASET(Layout_Cell) form.
  *  @see               PBblas/Types.Layout_Cell
  */
  EXPORT polyKernel(DATASET(PBblas.Types.Layout_cell) X, REAL deg=1.0) := FUNCTION
    value_t := PBblas.Types.value_t;
    dimension_t := PBblas.Types.dimension_t;
    value_t addOne(value_t v, dimension_t r, dimension_t c) := v + 1;
    value_t powerIt(value_t v, dimension_t r, dimension_t c) := POWER(v,deg);
    XXt := PBblas.gemm(FALSE, TRUE, 1.0, X, X);	
    XXtp1 := PBblas.Apply2Elements(XXt, addOne);
    XXtp1powDeg := PBblas.Apply2Elements(XXtp1, powerIt);
    RETURN XXtp1powDeg;
END;