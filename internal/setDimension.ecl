/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;

strict := FALSE;

EXPORT setDimension(DATASET(PBblas.Types.Layout_cell) A, UNSIGNED I, UNSIGNED J) := 
 IF(strict,
    IF(EXISTS(A(x=I,y=J)), A(x<=I,y<=J), 
    A(x<=I,y<=J)+DATASET([{1,I,J,0}], PBblas.Types.Layout_cell)),A);