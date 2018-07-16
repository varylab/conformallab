ref = 42.04533007815184

clusteringRandom <- abs(c(42.50382781804303,41.974235315265574) - ref)
clusteringRandom_range <- c(1.224126812336549,0.8608041541590538)
clusteringFibonacci <- abs(c(41.9191707121,42.0781627382468) - ref)
clusteringFibonacci_range <- abs(c(1.3459407174175415,0.9132731645532122) - ref)
homogeneousRandom <- abs(c(43.082974780281376,41.95100308725416) - ref)
homogeneousRandom_range <- abs(c(1.1930099638054783,0.5916161894882757) - ref)
homogeneousFibonacci <- abs(c(43.013935360062035,41.957810269492946) - ref)
homogeneousFibonacci_range <- abs(c(1.2305330240830745,0.5622312504095085) - ref)

plot(y=clusteringRandom, ylab="Error", x=clusteringRandom_range, xlab="Max Edge Length", log="xy", type="o", pch=22, lty=1, col=1
lines(y=clusteringFibonacci, x=clusteringFibonacci_range, type="o", pch=22, lty=1, col=2, ylim=g_range)
lines(y=homogeneousRandom, x=homogeneousRandom_range, type="o", pch=22, lty=1, col=3, ylim=g_range)
lines(y=homogeneousFibonacci, x=homogeneousFibonacci_range, type="o", pch=22, lty=1, col=4, ylim=g_range)

legend("topright", c("Clustering Random","Clustering Fibonacci","Homogeneous Random","Homogeneous Fibonacci"), col=c(1,2,3,4), lty=1);
title("General Torus")

# ---------------------------
#  REF  : ((-6.390807776698908+1.0967707323942397i))
# |REF| : 42.04533007815184
# Clustering Random ---------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.408161390354468+1.1997063862517867i))
# |P[0]|: 42.50382781804303
# |L[0]|: 1.224126812336549
# |a[0]|: 0.73909587319825
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.387137825839307+1.085682139946265i))
# |P[1]|: 41.974235315265574
# |L[1]|: 0.8608041541590538
# |a[1]|: 0.16223931793765595
# Clustering Fibonacci ------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.377920022886106+1.1141395306551558i))
# |P[0]|: 41.9191707121
# |L[0]|: 1.3459407174175415
# |a[0]|: 0.28799954735074484
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.392963599876307+1.0991720288032167i))
# |P[1]|: 42.0781627382468
# |L[1]|: 0.9132731645532122
# |a[1]|: 0.29551360991985454
# Homogeneous Random --------
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.4612103466119954+1.15574029830039i))
# |P[0]|: 43.082974780281376
# |L[0]|: 1.1930099638054783
# |a[0]|: 0.6707220103910131
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.376006597727362+1.139097429103986i))
# |P[1]|: 41.95100308725416
# |L[1]|: 0.5916161894882757
# |a[1]|: 0.49447034183902866
# Homogeneous Fibonacci -----
#  S[0] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 84, oriented edges: 504, faces: 168]
#  P[0] : ((-6.466739411433853+1.0932596876635423i))
# |P[0]|: 43.013935360062035
# |L[0]|: 1.2305330240830745
# |a[0]|: 0.8310567008628671
#  S[1] : CoHDS<CoVertex,CoEdge,CoFace>[vertices: 804, oriented edges: 4824, faces: 1608]
#  P[1] : ((-6.386325758267598+1.0828913051551996i))
# |P[1]|: 41.957810269492946
# |L[1]|: 0.5622312504095085
# |a[1]|: 0.529482945935935
