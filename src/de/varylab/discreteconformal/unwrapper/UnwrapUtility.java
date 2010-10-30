package de.varylab.discreteconformal.unwrapper;

import static de.varylab.discreteconformal.unwrapper.UnwrapUtility.BoundaryMode.Isometric;
import static de.varylab.discreteconformal.unwrapper.UnwrapUtility.QuantizationMode.AllAngles;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.signum;
import geom3d.Basis;
import geom3d.Point;
import geom3d.Vector;

import java.util.Collection;

import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.Adapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureField;
import de.jtem.halfedgetools.adapter.type.Normal;
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

	
	public static enum QuantizationMode {
		AllAngles,
		Quads,
		QuadsStrict,
		Hexagons,
		Straight
	}
	
	
	public static enum BoundaryMode {
		Isometric,
		Conformal,
		ConformalCurvature
	}
	
	
	
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
		for (final CoEdge e : hds.getPositiveEdges()) {
			final double l = e.getLength();
			e.setLambda(log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		Collection<CoVertex> vertexBoundary = HalfEdgeUtils.boundaryVertices(hds);
		Collection<CoEdge> edgeBoundary = HalfEdgeUtils.boundaryEdges(hds);
		
		// find reference normal
		Vector refNormal = new Vector(0, 0, 0);
		Point middle = new Point();
		for (CoVertex v : vertexBoundary) {
			middle.add(v.getPosition());
		}
		middle.times(1.0 / vertexBoundary.size());
		for (CoEdge e : edgeBoundary) {
			Point p1 = e.getStartVertex().getPosition();
			Point p2 = e.getTargetVertex().getPosition();
			Vector bv = p1.vectorTo(p2);
			Vector mv = p1.vectorTo(middle);
			Vector n = new Vector(mv).cross(bv).normalize();
			refNormal.add(n);
		}
		refNormal.times(1.0 / edgeBoundary.size());
		refNormal.normalize();
		
		// set thetas and solver indices
		int bSize = 0;
		int dim = 0;
		double bSum = 0.0;
		for (CoVertex v : hds.getVertices()) {
			if (!HalfEdgeUtils.isBoundaryVertex(v)) {
				if (v.getCustomInfo() != null && v.getCustomInfo().useCustomTheta) {
					bSum += 2 * PI - v.getCustomInfo().theta;
					v.setTheta(v.getCustomInfo().theta);
				} else {
					v.setTheta(2 * PI);
				}
				v.setSolverIndex(dim++);
				continue;
			}
			BoundaryMode mode = bm;
			if (v.getCustomInfo() != null) {
				if (v.getCustomInfo().useCustomTheta) {
					mode = BoundaryMode.Conformal;
				} else {
					mode = v.getCustomInfo().boundaryMode;
				}
			}
			switch (mode) {
			case Isometric:
				v.setTheta(0.0);
				v.setSolverIndex(-1);
				break;
			case Conformal:
				double theta = 0;
				if (v.getCustomInfo() != null) {
					if (v.getCustomInfo().useCustomTheta) {
						v.setTheta(v.getCustomInfo().theta);
						v.setSolverIndex(dim++);
						bSize++;
						bSum += Math.PI - v.getTheta();
						continue;
					}
				}
				for (CoEdge edge : HalfEdgeUtils.incomingEdges(v)) {
					if (edge.getLeftFace() == null) {
						Point pv = v.getPosition();
						Point p1 = edge.getStartVertex().getPosition();
						Point p2 = edge.getNextEdge().getTargetVertex().getPosition();
						Vector v1 = pv.vectorTo(p1);
						Vector v2 = pv.vectorTo(p2);
						Vector cr = new Vector(v1).cross(v2);
						double x = v1.dot(v2);
						double y = cr.getLength();
						theta = Math.atan2(y, x);
						double check = cr.dot(refNormal);
						if (check < 0) {
							theta = 2*PI - theta;
						}
						QuantizationMode qMode = qm;
						if (v.getCustomInfo() != null) {
							qMode = v.getCustomInfo().quantizationMode;
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
			CoFace f0 = hds.getFace(0);
			Vector f01 = f0.getBoundaryEdge().getVector();
			Vector f02 = f0.getBoundaryEdge().getPreviousEdge().getOppositeEdge().getVector();
			Vector f0N = new Vector(aSet.get(Normal.class, f0, double[].class));
			Basis f0B = new Basis(f01, f02, f0N);
			double surfaceOrientation = signum(f0B.getDeterminant());
			System.out.println("surface orientation is " + surfaceOrientation);
			Adapter<double[]> cVec = aSet.query(CurvatureField.class, hds.getEdgeClass(), double[].class);
			if (cVec == null) throw new RuntimeException("No curvature vector field on edges found");
			Collection<CoEdge> bList = HalfEdgeUtils.boundaryEdges(hds);
			CoEdge e = bList.iterator().next();
			CoEdge e0 = e;
			
			Double alphaP = null;
			Double firstAlpha = null;
			do {
				CoVertex v = e.getStartVertex();
				Vector eVec = e.getVector();
				double[] xArr = cVec.get(e, aSet);
				double[] xOppArr = cVec.get(e.getOppositeEdge(), aSet);
				Vector xVec = new Vector(xArr);
				double[] nArr = aSet.get(Normal.class, e.getRightFace(), double[].class);
				Vector nVec = new Vector(nArr);
				Basis B = new Basis(nVec, eVec, xVec);
				double sign = Math.signum(B.getDeterminant());
				double alpha = normalizeAngle(eVec.getAngle(xVec) * sign);
				
				if (alphaP != null) { // check if we must flip viVec
					double th = 0.0;
					for (CoEdge edge : HalfEdgeUtils.incomingEdges(v)) {
						if (edge.getLeftFace() == null) continue;
						th += getAngle(edge);
					}
					if (th < PI) {
						System.out.println("corner");
					}
					th *= surfaceOrientation;
					double gamma =  normalizeAngle(th - alphaP + alpha);
					if(abs(gamma) > PI/2) { // flip x
						System.out.println("flip at " + v.getIndex());
						alpha += PI;
						normalizeAngle(alpha);
						Rn.times(xArr, -1, xArr); // flip vector
						Rn.times(xOppArr, -1, xOppArr); // flip opposite vector
					}
					double a1 = alpha < 0 ? alpha + 2*PI : alpha;
					double a2 = alphaP < 0 ? alphaP + 2*PI : alphaP;
					double theta = a1 - a2;
					if (theta > 0) {
						theta = 2*PI - theta;
					}
					theta = abs(theta);
					System.out.println("sum theta: " + th);
					System.out.println("gamma: " + gamma);
					System.out.println("a: " + alpha + "   ap: " + alphaP);
					System.out.println("Theta at " + v.getIndex() + ": " + theta);
					System.out.println("");
					bSum += Math.PI - theta;
					v.setTheta(theta);
					v.setSolverIndex(dim++);
					if (e == e0) break;
				} else {
					firstAlpha = alpha;
				}
				if (alpha <= 0) {
					alphaP = alpha + PI;
				} else {
					alphaP = alpha - PI;
				}
				e = e.getNextEdge();
			} while (true);
			if (firstAlpha != alphaP) {
				throw new RuntimeException("Parallel transport along the boundary was discontinous!");
			}
		}
		
		// make conformal boundary embeddable
		if (bm == BoundaryMode.Conformal || bm == BoundaryMode.ConformalCurvature) {
//			double factor = (bSize - 2) * PI / bSum;
			System.out.println("Gauss-Bonnet sum: " + (bSum / PI));
//			for (CoVertex v : vertexBoundary) {
//				double theta = v.getTheta();
//				v.setTheta(theta * factor);
//			}
		}
		
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
	public static double getAngle(CoEdge e) {
		CoVertex v = e.getTargetVertex();
		Point pv = v.getPosition();
		Point p1 = e.getStartVertex().getPosition();
		Point p2 = e.getNextEdge().getTargetVertex().getPosition();
		Vector v1 = pv.vectorTo(p1);
		Vector v2 = pv.vectorTo(p2);
		Vector cr = new Vector(v1).cross(v2);
		double x = v1.dot(v2);
		double y = cr.getLength();
		return Math.atan2(y, x);
	}
	
	
	
	
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public static int prepareInvariantDataHyperbolic(CoHDS hds) {
		// set initial lambdas
		for (final CoEdge e : hds.getPositiveEdges()) {
			final double l = e.getLength();
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
