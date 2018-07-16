ref = 42.045330078150876

clusteringRandom <- abs(c(42.197623206842685,41.95666527624924,42.03399370173546) - ref)
clusteringRandom_range <- c(1.3079826509274988,0.8355187135640042,1.0010223212024494)
clusteringFibonacci <- abs(c(42.34430232958934,42.078162738242796,42.04819445017103) - ref)
clusteringFibonacci_range <- c(1.345940717417541,1.0109661099940133,1.0223008048195334)
homogeneousRandom <- abs(c(42.499454302805525,42.137532085731834,42.13904435963908) - ref)
homogeneousRandom_range <- c(1.124501342509803,0.5629643069546785,0.18588440987867938)
homogeneousFibonacci <- abs(c(42.19085233977455,42.317298113038895,42.172426046418124) - ref)
homogeneousFibonacci_range <- c(1.1463963718381056,0.49145895673975853,0.19653979825798895)

xlim = range(clusteringRandom_range, clusteringFibonacci_range, homogeneousRandom_range, homogeneousFibonacci_range)
ylim = range(clusteringRandom, clusteringFibonacci, homogeneousRandom, homogeneousFibonacci)

plot(
	y=clusteringRandom, ylab="Error", ylim=ylim,
	x=clusteringRandom_range, xlab="Max Edge Length", xlim=xlim,
	log="xy", type="o", pch=22, lty=1, col=1
)
lines(y=clusteringFibonacci, x=clusteringFibonacci_range, type="o", pch=22, lty=1, col=2)
lines(y=homogeneousRandom, x=homogeneousRandom_range, type="o", pch=22, lty=1, col=3)
lines(y=homogeneousFibonacci, x=homogeneousFibonacci_range, type="o", pch=22, lty=1, col=4)

legend("topright", c("Clustering Random","Clustering Fibonacci","Homogeneous Random","Homogeneous Fibonacci"), col=c(1,2,3,4), lty=1);
title("General Torus")# ---------------------------
#  REF  : ((-6.390807776698672+1.0967707323951763i))
# |REF| : 42.045330078150876
# Clustering Random ---------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.3981424716640865+1.1231189247496702i))
# |P[0]|: 42.197623206842685
# |L[0]|: 1.3079826509274988
# |a[0]|: 0.8089683400026558
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.384917388386183+1.0906398212668549i))
# |P[1]|: 41.95666527624924
# |L[1]|: 0.8355187135640042
# |a[1]|: 0.19560922932189378
#  S[2] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 8004, oriented edges: 48024, faces: 16008]
#  P[2] : ((-6.389998973949392+1.0963151073761515i))
# |P[2]|: 42.03399370173546
# |L[2]|: 1.0010223212024494
# |a[2]|: 0.09684593715702795
# Clustering Fibonacci ------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.403038495687914+1.1599139417767128i))
# |P[0]|: 42.34430232958934
# |L[0]|: 1.345940717417541
# |a[0]|: 0.6539187970619849
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.392963599876142+1.0991720288023512i))
# |P[1]|: 42.078162738242796
# |L[1]|: 1.0109661099940133
# |a[1]|: 0.23854151009830143
#  S[2] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 8004, oriented edges: 48024, faces: 16008]
#  P[2] : ((-6.390926558324666+1.0973842427661682i))
# |P[2]|: 42.04819445017103
# |L[2]|: 1.0223008048195334
# |a[2]|: 0.0813020699348767
# Homogeneous Random --------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.430741263806135+1.0700565876562749i))
# |P[0]|: 42.499454302805525
# |L[0]|: 1.124501342509803
# |a[0]|: 0.5185189272918017
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.4000388154749+1.084912552303855i))
# |P[1]|: 42.137532085731834
# |L[1]|: 0.5629643069546785
# |a[1]|: 0.4243758287590042
#  S[2] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 8004, oriented edges: 48024, faces: 16008]
#  P[2] : ((-6.397671768399933+1.099472738842376i))
# |P[2]|: 42.13904435963908
# |L[2]|: 0.18588440987867938
# |a[2]|: 0.30342472852595054
# Homogeneous Fibonacci -----
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.395881860921408+1.1329375803242907i))
# |P[0]|: 42.19085233977455
# |L[0]|: 1.1463963718381056
# |a[0]|: 0.6505368703305543
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.409679661792427+1.110542545895211i))
# |P[1]|: 42.317298113038895
# |L[1]|: 0.49145895673975853
# |a[1]|: 0.49447034183906247
#  S[2] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 8004, oriented edges: 48024, faces: 16008]
#  P[2] : ((-6.400761477468126+1.0966670210132994i))
# |P[2]|: 42.172426046418124
# |L[2]|: 0.19653979825798895
# |a[2]|: 0.24139553390661117
