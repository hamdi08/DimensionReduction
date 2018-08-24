IMPORT PBblas;
IMPORT $.^ as dr;
IMPORT $ as Tests;
IMPORT Tests.DiffReport as diffRep;


A := DATASET([
              {1, 1, 1, -1.0}, {1, 1, 2, -1.0}, {1, 1, 3, 1.0},
              {1, 2, 1, 1.0}, {1, 2, 2, 3.0}, {1, 2, 3, 3.0},
              {1, 3, 1, -1.0}, {1, 3, 2, -1.0}, {1, 3, 3, 5.0},
              {1, 4, 1, 1.0}, {1, 4, 2, 3.0}, {1, 4, 3, 7.0}],
              PBblas.Types.Layout_Cell);
QR := dr.QR(A);
qr_comp_id := ENUM (Q = 1,  R = 2);
Q := dr.internal.PBMU.From(QR, qr_comp_id.Q);
R := dr.internal.PBMU.From(QR, qr_comp_id.R);

A_qr := PBblas.gemm(FALSE, FALSE, 1.0, Q, R);
//OUTPUT(Q);
//OUTPUT(R);
test := diffRep.Compare_Cells('TEST -- QR decom -- A = QR', A, A_qr);
OUTPUT(test);
//EXPORT qrTest := test;