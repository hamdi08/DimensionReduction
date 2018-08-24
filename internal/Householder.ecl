/*##############################################################################
## HPCC SYSTEMS software Copyright (C) 2018 HPCC Systems.  All rights reserved.
############################################################################## */
IMPORT PBblas;
/**
  *  This module does the Householder transformation. 
  *  <p> Computes the Householder reflection matrix for use with QR decomposition.
  
  *  @param vec    Vector or column matrix.
  *  @param k      index < length(vec).
  *  @return H     matrix H that annihilates entries in the product H*vec below index k.
  *  @see          PBblas/Types.Layout_Cell
  */
EXPORT Householder(DATASET(PBblas.Types.Layout_Cell) vec, UNSIGNED k, UNSIGNED Dim=1) := MODULE
 EXPORT HTA := MODULE
  EXPORT Default := MODULE, VIRTUAL
   EXPORT IdentityM := IF(Dim > COUNT(vec), PBblas.eye(Dim), PBblas.eye(COUNT(vec))); 
   //COUNT(X) --> length of the vector
   EXPORT DATASET(PBblas.Types.Layout_Cell) hReflection := DATASET([], PBblas.Types.Layout_Cell);
  END;
  // Householder Vector
  HouseV(DATASET(PBblas.Types.Layout_Cell) vec, UNSIGNED k) := FUNCTION
   xk := vec(x=k)[1].v;
   alpha := IF(xk >= 0, -1 * PBblas.vec.Norm(vec), PBblas.vec.Norm(vec));
   vk := IF (alpha=0, 1, SQRT(0.5*(1-xk/alpha)));
   p := - alpha * vk;
   RETURN PROJECT(vec, TRANSFORM(PBblas.Types.Layout_Cell, SELF.v := IF(LEFT.x=k, vk, LEFT.v/(2*p)), SELF :=LEFT));
  END;
  // Source: Atkinson, Section 9.3, p. 611	
  EXPORT Atkinson := MODULE(Default)
   hV := HouseV(vec(x>=k),k);
   houseVec := PBblas.vec.ToCol(hV, 1);
   EXPORT DATASET(PBblas.Types.Layout_Cell) hReflection := PBblas.axpy(-2, PBblas.gemm(FALSE, TRUE, 1.0, houseVec, houseVec), IdentityM);
  END;
  // Source: Golub and Van Loan, "Matrix Computations" p. 210
  EXPORT Golub := MODULE(Default)
   VkValue := vec(x=k)[1].v;
   VkPlus := vec(x>k);
   sigma := PBblas.vec.Dot(VkPlus, VkPlus);
   mu := SQRT(VkValue*VkValue + sigma);
   newVkValue := IF(sigma=0,1,IF(VkValue<=0, VkValue-mu, -1 * sigma/(VkValue+mu) ));
   beta := IF( sigma=0, 0, 2*(newVkValue*newVkValue)/(sigma + (newVkValue*newVkValue)));
   newVkElem0 := vec[1];
   newVkElem := PROJECT(newVkElem0,TRANSFORM(PBblas.Types.Layout_cell, SELF.x:=k, SELF.y:=1,SELF.v := newVkValue));
   hV := PROJECT(newVkElem + VkPlus,TRANSFORM(PBblas.Types.Layout_cell,SELF.v:=LEFT.v/newVkValue, SELF := LEFT));
   EXPORT DATASET(PBblas.Types.Layout_cell) hReflection := PBblas.axpy(-1 * Beta, PBblas.gemm(FALSE, TRUE, 1.0, hV, hV), IdentityM);
  END;
 END;
 EXPORT Reflection(HTA.Default Control = HTA.Golub) := FUNCTION
  RETURN Control.hReflection;
 END;
END;