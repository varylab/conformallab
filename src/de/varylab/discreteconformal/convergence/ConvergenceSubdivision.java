package de.varylab.discreteconformal.convergence;

import java.io.FileWriter;
import java.util.HashSet;
import java.util.Set;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.algorithm.subdivision.LoopLinear;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.calculator.SubdivisionCalculator;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class ConvergenceSubdivision extends ConvergenceSeries {

	private int 
		maxSubdivision = 0,
		numExtraPoints = 0;
	private LoopLinear 
		loop = new LoopLinear();
	private SubdivisionCalculator 
		sc = new SubdivisionCalculator();
	
	
	public ConvergenceSubdivision() {
	}

	public ConvergenceSubdivision(
		double[][] vertices, 
		int[] branchIndices,
		Complex tauExpected, 
		FileWriter errorWriter, 
		int maxSubdivision,
		LoopLinear loop, 
		SubdivisionCalculator sc
	) {
		super(vertices, branchIndices, tauExpected, errorWriter);
		this.maxSubdivision = maxSubdivision;
		this.loop = loop;
		this.sc = sc;
	}


	@Override
	protected OptionSet configureAndParseOptions(OptionParser p, String... args) {
		OptionSpec<Integer> mixSubdivisionSpec = p.accepts("max", "Maxmum number of subdivision steps").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> numExtraPointsSpec = p.accepts("extra", "Number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSet result = p.parse(args);
		maxSubdivision = mixSubdivisionSpec.value(result);
		numExtraPoints = numExtraPointsSpec.value(result);
		return result;
	}
	
	@Override
	protected void perform() throws Exception {
		writeComment("numVertex[1], absErr[2], argErr[3], reErr[4], imErr[5], gradNormSq[6]");
		for (int i = 0; i < maxSubdivision; i ++) {
			CoHDS hds = new CoHDS();
			// predefined vertices
			for (int vi = 0; vi < vertices.length; vi++) {
				CoVertex v = hds.addNewVertex();
				v.getPosition().set(vertices[vi][0], vertices[vi][1], vertices[vi][2]);	
			}
			// extra points
			for (int j = 0; j < numExtraPoints; j++) {
				CoVertex v = hds.addNewVertex();
				v.getPosition().set(rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian());
				v.getPosition().normalize();
			}
			ConvexHull.convexHull(hds, sc, 1E-8);
			// subdivision
			for (int si = 0; si < i; si++) {
				CoHDS subdivided = new CoHDS();
				loop.subdivide(hds, subdivided, sc, sc, sc);
				hds = subdivided;
				// project to the sphere in every step
				for (CoVertex v : hds.getVertices()) {
					v.getPosition().normalize();
				}
			}
			int numVerts = hds.numVertices();
			Complex tau = null;
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				DiscreteEllipticUtility.generateEllipticImage(hds, 0, glueSet, branchIndices);
				tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(hds, 1E-8);
			} catch (Exception e) {
				System.out.println("Error: " + e.getMessage());
				continue;
			}
			double absErr = tau.abs() - getExpectedTau().abs();
			double argErr = tau.arg() - getExpectedTau().arg();
			double reErr = tau.re - getExpectedTau().re;
			double imErr = tau.im - getExpectedTau().im;
			writeErrorLine(numVerts + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + EuclideanUnwrapperPETSc.lastGNorm);
		}
	}

}
