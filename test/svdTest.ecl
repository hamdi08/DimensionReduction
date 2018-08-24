IMPORT PBblas;
IMPORT $.^ as dr;
IMPORT $ as Tests;
IMPORT Tests.DiffReport as diffRep;
//Symmetric matrix
A := DATASET([
              {1, 1, 1, 1.0}, {1, 1, 2, 0.7355303763393315},
              {1, 2, 1, 0.7355303763393315}, {1, 2, 2, 1.0}],PBblas.Types.Layout_Cell);
U := dr.svd(A).UComp;
S := dr.svd(A).SComp;
V := dr.svd(A).VComp;

US := PBblas.gemm(FALSE, FALSE, 1.0, U, S);
USV_t := PBblas.gemm(FALSE, TRUE, 1.0, US, V);

test := diffRep.Compare_Cells('TEST -- SVD -- A = USV_t', A, USV_t);
EXPORT svdTest := test;