IMPORT PBblas;
IMPORT $.^ as dr;
IMPORT $ as Tests;
IMPORT Tests.DiffReport as diffRep;
//A is a matrix where both columns are maximally redundant
A := DATASET([
              {1, 1, 1, 1.0}, {1, 1, 2, 2.0},
              {1, 2, 1, 3.0}, {1, 2, 2, 4.0},
              {1, 3, 1, 5.0}, {1, 3, 2, 6.0},
              {1, 4, 1, 7.0}, {1, 4, 2, 8.0}],
              PBblas.Types.Layout_Cell);
cov := dr.cov(A);
//OUTPUT(cov, NAMED('cov'));
A_dims := PBblas.internal.matdims.FromCells(A);
n := A_dims[1].m_cols;
onesN := dr.internal.ones(n,n);
//OUTPUT(onesN, NAMED('oneN'));
test := diffRep.Compare_Cells('TEST -- Covmat --  (1/nRows)*A_t*A', onesN, cov);
//OUTPUT(test);
EXPORT covTest := test;
