ref = 42.04533007815184

clusteringRandom <- abs(c(42.50382781804303,41.97423531528474) - ref)
clusteringRandom_range <- c(1.224126812336549,0.9152847160679753)
clusteringFibonacci <- abs(c(42.81214431868112,42.01516187276862) - ref)
clusteringFibonacci_range <- c(1.3459407174175415,0.8900560669847543)
homogeneousRandom <- abs(c(42.20530460196575,42.08525096129623) - ref)
homogeneousRandom_range <- c(1.1930099638054783,0.5916161894882757)
homogeneousFibonacci <- abs(c(41.50758383039661,41.760763421760224) - ref)
homogeneousFibonacci_range <- c(1.2305330240830745,0.5622312504095085)

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
#  REF  : ((-6.390807776698908+1.0967707323942397i))
# |REF| : 42.04533007815184
# Clustering Random ---------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.408161390354468+1.1997063862517867i))
# |P[0]|: 42.50382781804303
# |L[0]|: 1.224126812336549
# |a[0]|: 0.73909587319825
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.3871378258406235+1.0856821399473464i))
# |P[1]|: 41.97423531528474
# |L[1]|: 0.9152847160679753
# |a[1]|: 0.19340347430477656
# Clustering Fibonacci ------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.44549414451179+1.1259438492861662i))
# |P[0]|: 42.81214431868112
# |L[0]|: 1.3459407174175415
# |a[0]|: 0.5019588503037915
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.388175925179275+1.0983488615774761i))
# |P[1]|: 42.01516187276862
# |L[1]|: 0.8900560669847543
# |a[1]|: 0.28315157722421724
# Homogeneous Random --------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.381113058915564+1.2193033795175847i))
# |P[0]|: 42.20530460196575
# |L[0]|: 1.1930099638054783
# |a[0]|: 0.6707220103910131
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.383924151347021+1.153587184024013i))
# |P[1]|: 42.08525096129623
# |L[1]|: 0.5916161894882757
# |a[1]|: 0.49447034183902866
# Homogeneous Fibonacci -----
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.352497435949849+1.0739460678484745i))
# |P[0]|: 41.50758383039661
# |L[0]|: 1.2305330240830745
# |a[0]|: 0.8310567008628671
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.371312030241188+1.0803455165196627i))
# |P[1]|: 41.760763421760224
# |L[1]|: 0.5622312504095085
# |a[1]|: 0.529482945935935
