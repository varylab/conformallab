package de.varylab.discreteconformal.util;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.signum;

import java.util.Collection;
import java.util.List;
import java.util.logging.Logger;

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
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.functional.ConformalFunctional;
import de.varylab.discreteconformal.functional.FunctionalAdapters.InitialEnergy;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin.TargetGeometry;
import de.varylab.discreteconformal.unwrapper.BoundaryMode;
import de.varylab.discreteconformal.unwrapper.QuantizationMode;
import de.varylab.discreteconformal.unwrapper.numerics.SimpleEnergy;

public class UnwrapUtility {

	private static Logger logger = Logger.getLogger(UnwrapUtility.class.getName());
	
	public static class ZeroInitialEnergy implements InitialEnergy<CoFace> {
		@Override
		public double getInitialEnergy(CoFace f) {
			return 0.0;
		}
	}
	
	public static class ZeroU implements DomainValue {
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
	public static int prepareInvariantDataEuclidean(ConformalFunctional<CoVertex, CoEdge, CoFace> fun, CoHDS hds, AdapterSet aSet) {
		return prepareInvariantDataEuclidean(fun, hds, BoundaryMode.Isometric, QuantizationMode.AllAngles, aSet);
	}
	
	
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public static int prepareInvariantDataEuclidean(ConformalFunctional<CoVertex, CoEdge, CoFace> fun, CoHDS hds, BoundaryMode bm, QuantizationMode qm, AdapterSet aSet) {
		// set initial lambdas
		TypedAdapterSet<double[]> da = aSet.querySet(double[].class);
		for (final CoEdge e : hds.getPositiveEdges()) {
			Double l = aSet.get(Length.class, e, Double.class);
			double lambda = fun.getLambda(l);
			e.setLambda(lambda);
			e.getOppositeEdge().setLambda(e.getLambda());
		}

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
						for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
							if (e.getLeftFace() == null) continue;
							theta += (getAngle(e, aSet) % PI);
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
							if((theta % (PI/2)) != 0) {
								logger.info(v + ": theta quantized to " + theta);
							}
							break;
						case QuadsStrict:
							q = theta % (PI/2);
							if (q > PI/4) {
								theta += (PI/2 - q);
							} else {
								theta -= q;
							}
							if((theta % PI) != 0) {
								logger.info(v + ": theta quantized to " + theta);
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
			default:
				break;
			}
		}
		
		if (bm == BoundaryMode.ConformalCurvature && bSize != 0) {
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
					double gamma =  normalizeAngle(th - alphaP + alpha);
					if(abs(gamma) > PI/2) { // flip x
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
		logger.info("Gauss-Bonnet sum: " + (bSum / PI) + " \u03C0");

		// circular quad edges
		for (CoEdge e : hds.getPositiveEdges()) {
			if (e.info != null && e.info.circularHoleEdge) {
				int index = dim++;
				e.setSolverIndex(index);
				e.getOppositeEdge().setSolverIndex(index);
				e.setPhi(e.info.phi);
				e.getOppositeEdge().setPhi(e.info.phi);
			} else {
				e.setSolverIndex(-1);
				e.getOppositeEdge().setSolverIndex(-1);
			}
		}
		
		// initial Euclidean energy
		ZeroU zeroU = new ZeroU();
		SimpleEnergy E = new SimpleEnergy();
		ZeroInitialEnergy zeroEnergy = new ZeroInitialEnergy();
		for (final CoFace f : hds.getFaces()) {
			E.setZero();
			fun.triangleEnergyAndAlphas(zeroU, f, E, zeroEnergy);
			f.setInitialEnergy(E.get());
		}
		return dim;
	}
	
	/**
	 * Normalize a given angle to the interval [-pi,pi]
	 * @param a
	 * @return
	 */
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
	public static int prepareInvariantDataHyperbolicAndSpherical(
		ConformalFunctional<CoVertex, CoEdge, CoFace> fun, 
		CoHDS hds, 
		AdapterSet a,
		DomainValue initialU
	) {
		return prepareInvariantDataHyperbolicAndSpherical(fun, hds, a, initialU, 1.0);
	}
	
	/**
	 * Compute algorithm invariant data
	 * @param boundary the boundary vertices which do not belong to the solver system
	 * @return the dimension of the parameter space
	 */
	public static int prepareInvariantDataHyperbolicAndSpherical(
		ConformalFunctional<CoVertex, CoEdge, CoFace> fun, 
		CoHDS hds, 
		AdapterSet a,
		DomainValue initialU,
		double scale
	) {
		// set initial lambdas
		for (final CoEdge e : hds.getPositiveEdges()) {
			try {
				Double l = a.get(Length.class, e, Double.class);
				l*= scale;
				double lambda = fun.getLambda(l);
				e.setLambda(lambda);
				e.getOppositeEdge().setLambda(lambda);
			} catch (Exception e1) {
				System.err.println("error setting lambdas in prepareInvariantData(): " + e1);
				return -1;
			}
		}
		
		// set thetas and solver indices
		int dim = 0;
		for (final CoVertex v : hds.getVertices()) {
			if (v.info != null && v.info.useCustomTheta) {
				v.setTheta(v.info.theta);
			} else {
				v.setTheta(2 * PI);
			}
			v.setSolverIndex(dim++);
		}
		
		// circular quad edges
		for (CoEdge e : hds.getPositiveEdges()) {
			if (e.info != null && e.info.circularHoleEdge) {
				int index = dim++;
				e.setSolverIndex(index);
				e.getOppositeEdge().setSolverIndex(index);
			} else {
				e.setSolverIndex(-1);
				e.getOppositeEdge().setSolverIndex(-1);
			}
		}
		
		// initial hyperbolic energy
		ZeroInitialEnergy initialEnergy = new ZeroInitialEnergy();
		SimpleEnergy E = new SimpleEnergy();
		for (final CoFace f : hds.getFaces()) {
			E.setZero();
			fun.triangleEnergyAndAlphas(initialU, f, E, initialEnergy);
			f.setInitialEnergy(E.get());
		}
		return dim;
	}



	public static TargetGeometry calculateTargetGeometry(CoHDS hds) {
		int g = HalfEdgeUtils.getGenus(hds);
		List<List<CoEdge>> bc = HalfEdgeUtils.boundaryComponents(hds);
		return calculateTargetGeometry(g, bc.size());
	}


	public static TargetGeometry calculateTargetGeometry(int g, int numBoundaryComponents) {
		switch (g) {
		case 0:
			if (numBoundaryComponents > 0) {
				return TargetGeometry.Euclidean;
			} else {
				return TargetGeometry.Spherical;
			}
		case 1:
			return TargetGeometry.Euclidean;
		default:
			return TargetGeometry.Hyperbolic;
		}
	}


	public static void logSolutionStatus(Tao optimizer, Logger log) {
		log.info("TAO: " + optimizer.getSolutionStatus().reason);
		log.info("TAO iterate: " + optimizer.getSolutionStatus().iterate);
		log.info("TAO f: " + optimizer.getSolutionStatus().f);
		log.info("TAO gnorm: " + optimizer.getSolutionStatus().gnorm);
		log.info("TAO xdiff: " + optimizer.getSolutionStatus().xdiff);
	}
	
}