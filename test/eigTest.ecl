IMPORT PBblas;
IMPORT $.^ as dr;
IMPORT $ as Tests;
IMPORT Tests.DiffReport as diffRep;

A := DATASET([
              {1, 1, 1, 2.0}, {1, 1, 2, 3.0},
              {1, 2, 1, 2.0}, {1, 2, 2, 1.0}],PBblas.Types.Layout_Cell);
U := dr.eig(A).vectors;
S := dr.eig(A).valuesM;

//OUTPUT(U, NAMED('U'));
//OUTPUT(S, NAMED('S'));

US := PBblas.gemm(FALSE, FALSE, 1.0, U, S);
USU_t := PBblas.gemm(FALSE, TRUE, 1.0, US, U);

test := diffRep.Compare_Cells('TEST -- Eigendecomposition -- A = USU_t', A, USU_t);
//OUTPUT(test);

EXPORT eigTest := test;