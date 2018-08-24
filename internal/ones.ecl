/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
/**
  *  Returns a matrix of all ones with given shape.
  *  @param nr          Number of rows in unsigned.
  *  @param nc          Number of columns in unsigned.
  *  @return            Matrix in DATASET(Layout_Cell) form.
  *  @see               PBblas/Types.Layout_Cell
  */
  
EXPORT ones(UNSIGNED nr, UNSIGNED nc) := FUNCTION
 numIter := nc - 1;
 loopbody(DATASET(PBblas.Types.Layout_cell) m, UNSIGNED i) := FUNCTION
  m1 := PROJECT(m(y=i), TRANSFORM(PBblas.Types.Layout_cell,
                                  SELF.y:= i+1, SELF := LEFT));
  RETURN m+m1;
 END;
 RETURN LOOP(PBblas.vec.From(nr, 1), COUNTER <= numIter, loopbody(ROWS(LEFT), COUNTER) );
END;