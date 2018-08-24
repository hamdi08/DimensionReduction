IMPORT PBblas;
IMPORT $.^ as dr;
IMPORT $ as Tests;
IMPORT Tests.DiffReport as diffRep;


A := DATASET([
              {1, 1, 1, 7.0}, {1, 1, 2, 4.0}, {1, 1, 3, 3.0},
              {1, 2, 1, 4.0}, {1, 2, 2, 1.0}, {1, 2, 3, 8.0},
              {1, 3, 1, 6.0}, {1, 3, 2, 3.0}, {1, 3, 3, 5.0},
              {1, 4, 1, 8.0}, {1, 4, 2, 6.0}, {1, 4, 3, 1.0},
              {1, 5, 1, 8.0}, {1, 5, 2, 5.0}, {1, 5, 3, 7.0},
              {1, 6, 1, 7.0}, {1, 6, 2, 2.0}, {1, 6, 3, 9.0},
              {1, 7, 1, 5.0}, {1, 7, 2, 3.0}, {1, 7, 3, 3.0},
              {1, 8, 1, 9.0}, {1, 8, 2, 5.0}, {1, 8, 3, 8.0},
              {1, 9, 1, 7.0}, {1, 9, 2, 4.0}, {1, 9, 3, 5.0},
              {1, 10, 1, 8.0}, {1, 10, 2, 2.0}, {1, 10, 3, 2.0}],
              PBblas.Types.Layout_Cell);
Agg := TABLE(A, {wi_id, y, var := VARIANCE(GROUP, v), mean := AVE(GROUP, v)}, wi_id, y);
A_norm := JOIN(A, Agg, 
                LEFT.wi_id=RIGHT.wi_id AND LEFT.y= RIGHT.y, 
                TRANSFORM(PBblas.Types.Layout_cell, 
                          SELF.v := (LEFT.v - RIGHT.mean)/SQRT(RIGHT.var), 
                          SELF := LEFT), 
                          LOOKUP);
d := 2;
Z := dr.pca(A,d).ZComp;
URed := dr.pca(A,d).UReduce;
//OUTPUT(Z, NAMED('Z'));
//OUTPUT(URed, NAMED('UReduce'));

A_approx := PBblas.gemm(FALSE, TRUE, 1.0, Z, URed);
//OUTPUT(A_approx, NAMED('A_approx'));
//OUTPUT(A_norm, NAMED('A_norm'));

test := diffRep.Compare_Cells('TEST -- PCA ', A_norm, A_approx);
//there will be reconstruction errors
//OUTPUT(test);
EXPORT pcaTest := test;

