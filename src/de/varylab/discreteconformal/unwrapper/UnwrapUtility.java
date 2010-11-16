package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.unwrapper.BoundaryMode.Isometric;
import static de.varylab.discreteconformal.unwrapper.QuantizationMode.AllAngles;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.signum;

import java.util.Collection;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.Adapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.TypedAdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureField;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.halfedgetools.functional.DomainValue;
import de.varylab.discreteconformal.functional.ConformalAdapters.InitialEnergy;
import de.varylab.discreteconformal.functional.ConformalEuclideanFunctional;
import de.varylab.discreteconformal.functional.ConformalHyperbolicFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.ConformalEnergy;

public class UnwrapUtility {

	
	private static class ZeroInitialEnergy implements InitialEnergy<CoFace> {
		@Override
		public double getInitialEnergy(CoFace f) {
			return 0.0;
		}
	}
	
	private static class ZeroU implements DomainValue {
		@Override
		public void add(int i, double value) {
		}
		@Override
		public void set(int i, double value) {
		}
		@Override
		public void setZero() {
		}
		@Override
		public double get(int i) {
			return 0.0;
		}
	}
	
	
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public static int prepareInvariantDataEuclidean(CoHDS hds, AdapterSet aSet) {
		return prepareInvariantDataEuclidean(hds, Isometric, AllAngles, aSet);
	}
	
	
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public static int prepareInvariantDataEuclidean(CoHDS hds, BoundaryMode bm, QuantizationMode qm, AdapterSet aSet) {
		// set initial lambdas
		TypedAdapterSet<double[]> da = aSet.querySet(double[].class);
		for (final CoEdge e : hds.getPositiveEdges()) {
			double l = aSet.get(Length.class, e, Double.class);
			e.setLambda(log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		Collection<CoVertex> vertexBoundary = HalfEdgeUtils.boundaryVertices(hds);
		Collection<CoEdge> edgeBoundary = HalfEdgeUtils.boundaryEdges(hds);
		
		// find reference normal
		double[] refNormal = new double[3];
		double[] middle = new double[4];
		for (CoVertex v : vertexBoundary) {
			Pn.dehomogenize(v.P, v.P);
			Rn.add(middle, v.P, middle);
		}
		Rn.times(middle, 1.0 / vertexBoundary.size(), middle);
		middle[3] = 1.0;
		for (CoEdge e : edgeBoundary) {
			double[] p = e.getStartVertex().P;
			Pn.dehomogenize(p, p);
			double[] bv = da.get(EdgeVector.class, e);
			double[] mv = {middle[0] - p[0], middle[1] - p[1], middle[2] - p[1]};
			double[] n = Rn.crossProduct(null, mv, bv);
			Rn.add(refNormal, n, refNormal);
		}
		Rn.normalize(refNormal, refNormal);
		
		// set thetas and solver indices
		int bSize = 0;
		int dim = 0;
		double bSum = 0.0;
		for (CoVertex v : hds.getVertices()) {
			if (!HalfEdgeUtils.isBoundaryVertex(v)) {
				if (v.info != null && v.info.useCustomTheta) {
					bSum += 2 * PI - v.info.theta;
					v.setTheta(v.info.theta);
				} else {
					v.setTheta(2 * PI);
				}
				v.setSolverIndex(dim++);
				continue;
			}
			BoundaryMode mode = bm;
			if (v.info != null) {
				if (v.info.useCustomTheta) {
					mode = BoundaryMode.Conformal;
				} else {
					mode = v.info.boundaryMode;
				}
			}
			switch (mode) {
			case Isometric:
				v.setTheta(0.0);
				v.setSolverIndex(-1);
				break;
			case Conformal:
				double theta = 0;
				if (v.info != null) {
					if (v.info.useCustomTheta) {
						v.setTheta(v.info.theta);
						v.setSolverIndex(dim++);
						bSize++;
						bSum += Math.PI - v.getTheta();
						continue;
					}
				}
				for (CoEdge edge : HalfEdgeUtils.incomingEdges(v)) {
					if (edge.getLeftFace() == null) {
						double[] v1 = da.get(EdgeVector.class, edge.getOppositeEdge());
						double[] v2 = da.get(EdgeVector.class, edge.getNextEdge());
						double[] cr = Rn.crossProduct(null, v1, v2); 
						double x = Rn.innerProduct(v1, v2);
						double y = Rn.euclideanNorm(cr);
						theta = Math.atan2(y, x);
						double check = Rn.innerProduct(cr, refNormal);
						if (check < 0) {
							theta = 2*PI - theta;
						}
						QuantizationMode qMode = qm;
						if (v.info != null) {
							qMode = v.info.quantizationMode;
						}
						switch (qMode) {
						case AllAngles: break;
						case Quads:
							double q = theta % (PI/4);
							if (q > PI/8) {
								theta += (PI/4 - q);
							} else {
								theta -= q;
							}
							break;
						case QuadsStrict:
							q = theta % (PI/2);
							if (q > PI/4) {
								theta += (PI/2 - q);
							} else {
								theta -= q;
							}
							break;							
						case Hexagons:
							q = theta % (PI/6);
							if (q > PI/12) {
								theta += (PI/6 - q);
							} else {
								theta -= q;
							}
							break;
						case Straight:
							theta = PI;
							break;
						}
					}
				}
				bSize++;
				bSum += Math.PI - theta;
				v.setTheta(theta);
				v.setSolverIndex(dim++);
				break;
			case ConformalCurvature:
				bSize++;
				break;
			}
		}
		
		if (bm == BoundaryMode.ConformalCurvature && bSize != 0) {
//			CoFace f0 = hds.getFace(0);
//			double[] f01 = da.get(EdgeVector.class, f0.getBoundaryEdge());
//			double[] f02 = da.get(EdgeVector.class, f0.getBoundaryEdge().getPreviousEdge().getOppositeEdge());
//			double[] f0N = da.get(Normal.class, f0);
//			double[][] f0B = {f01, f02, f0N};
//			double surfaceOrientation = signum(Rn.determinant(f0B));
//			System.out.println("surface orientation is " + surfaceOrientation);
			Adapter<double[]> cVec = da.query(CurvatureField.class, hds.getEdgeClass());
			if (cVec == null) throw new RuntimeException("No curvature vector field on edges found");
			Collection<CoEdge> bList = HalfEdgeUtils.boundaryEdges(hds);
			CoEdge e = bList.iterator().next();
			CoEdge e0 = e;
			
			Double alphaP = null;
			do {
				CoVertex v = e.getStartVertex();
				double[] eVec = da.get(EdgeVector.class, e);
				double[] xVec = cVec.get(e, aSet);
				double[] nVec = da.get(Normal.class, e.getRightFace());
				double[][] B = {nVec, eVec, xVec};
				double sign = Math.signum(Rn.determinant(B));
				double alpha = normalizeAngle(Rn.euclideanAngle(eVec, xVec) * sign);
				
				if (alphaP != null) { // check if we must flip viVec
					double th = 0.0;
					for (CoEdge edge : HalfEdgeUtils.incomingEdges(v)) {
						if (edge.getLeftFace() == null) continue;
						th += getAngle(edge, aSet);
					}
//					th *= surfaceOrientation;
					double gamma =  normalizeAngle(th - alphaP + alpha);
					if(abs(gamma) > PI/2) { // flip x
//						System.out.println("flip at " + v.getIndex());
						alpha += PI;
						normalizeAngle(alpha);
						Rn.times(xVec, -1, xVec); // flip vector
					}
					double a1 = alpha < 0 ? alpha + 2*PI : alpha;
					double a2 = alphaP < 0 ? alphaP + 2*PI : alphaP;
					double theta = a1 - a2;
					if (theta > 0) {
						theta = 2*PI - theta;
					}
					theta = abs(theta);
//					System.out.println("sum theta: " + th);
//					System.out.println("gamma: " + gamma);
//					System.out.println("a: " + alpha + "   ap: " + alphaP);
//					System.out.println("Theta at " + v.getIndex() + ": " + theta);
//					System.out.println("");
					bSum += Math.PI - theta;
					v.setTheta(theta);
					v.setSolverIndex(dim++);
					if (e == e0) {
						break;
					}
				}
				if (alpha <= 0) {
					alphaP = alpha + PI;
				} else {
					alphaP = alpha - PI;
				}
				e = e.getNextEdge();
			} while (true);
		}
		
		System.out.println("Gauss-Bonnet sum: " + (bSum / PI));
		
		// initial Euclidean energy
		ZeroU zeroU = new ZeroU();
		CVariable var = new CVariable();
		CLambda lambda = new CLambda();
		CAlpha alpha = new CAlpha();
		CTheta theta = new CTheta();
		ZeroInitialEnergy zeroEnergy = new ZeroInitialEnergy();
		ConformalEnergy E = new ConformalEnergy();
		ConformalEuclideanFunctional<CoVertex, CoEdge, CoFace> func = new ConformalEuclideanFunctional<CoVertex, CoEdge, CoFace>(var, theta, lambda, alpha, zeroEnergy);
		for (final CoFace f : hds.getFaces()) {
			E.setZero();
			func.triangleEnergyAndAlphas(hds, zeroU, f, E);
			f.setInitialEnergy(E.get());
		}
		return dim;
	}
	
	
	
	private static double normalizeAngle(double a) {
		a %= 2*PI;
		if (abs(a) > PI) {
			a -= signum(a) * 2*PI;
		}
		return a;
	}
	
	
	
	
	/**
	 * Calculate the angle at the target vertex of e
	 * @param e
	 * @return
	 */
	public static double getAngle(CoEdge e, AdapterSet a) {
		double[] v1 = a.get(EdgeVector.class, e.getOppositeEdge(), double[].class);
		double[] v2 = a.get(EdgeVector.class, e.getNextEdge(), double[].class);
		double[] cr = Rn.crossProduct(null, v1, v2);
		double x = Rn.innerProduct(v1, v2);
		double y = Rn.euclideanNorm(cr);
		return Math.atan2(y, x);
	}
	
	
	
	
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public static int prepareInvariantDataHyperbolic(CoHDS hds, AdapterSet a) {
		// set initial lambdas
		for (final CoEdge e : hds.getPositiveEdges()) {
			double l = a.get(Length.class, e, Double.class);
			e.setLambda(2*log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		
		// set thetas and solver indices
		int dim = 0;
		for (final CoVertex v : hds.getVertices()) {
			v.setTheta(2 * PI);
			v.setSolverIndex(dim++);
		}
		// initial hyperbolic energy
		ZeroU zeroU = new ZeroU();
		CVariable var = new CVariable();
		CLambda lambda = new CLambda();
		CTheta theta = new CTheta();
		CAlpha alpha = new CAlpha();
		ZeroInitialEnergy zeroEnergy = new ZeroInitialEnergy();
		ConformalEnergy E = new ConformalEnergy();
		ConformalHyperbolicFunctional<CoVertex, CoEdge, CoFace> func = new ConformalHyperbolicFunctional<CoVertex, CoEdge, CoFace>(var, theta, lambda, alpha, zeroEnergy);
		for (final CoFace f : hds.getFaces()) {
			E.setZero();
			func.triangleEnergyAndAlphas(zeroU, f, E);
			f.setInitialEnergy(E.get());
		}
		return dim;
	}
	
	
	
	
}
