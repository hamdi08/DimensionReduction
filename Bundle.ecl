IMPORT Std;
EXPORT Bundle := MODULE(Std.BundleBase)
 EXPORT Name := 'DimensionReduction';
 EXPORT Description := 'Dimensionality Reduction using PCA and kernel PCA';
 EXPORT Authors := ['HPCCSystems'];
 EXPORT License := 'http://www.apache.org/licenses/LICENSE-2.0';
 EXPORT Copyright := 'Copyright (C) 2018, 2017 HPCC Systems';
 EXPORT DependsOn := ['PBblas'];
 EXPORT Version := '2.0';
 EXPORT PlatformVersion := '6.2.0';
END;