package de.varylab.discreteconformal.convergence;

import java.io.FileOutputStream;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.math.Pn;
import de.jtem.halfedgetools.adapter.TypedAdapterSet;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.algorithm.subdivision.Loop;
import de.jtem.mfc.field.Complex;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.ConformalDataList;
import de.varylab.conformallab.data.types.HalfedgeEmbedding;
import de.varylab.conformallab.data.types.HalfedgeMap;
import de.varylab.conformallab.data.types.ObjectFactory;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class ConvergenceSubdivision extends ConvergenceSeries {

	private static Logger
		log = Logger.getLogger(ConvergenceSubdivision.class.getName());
	private int 
		maxSubdivision = 0,
		numExtraPoints = 0;
	private Loop 
		loop = new Loop();
	
	
	public ConvergenceSubdivision() {
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
		writeComment("numVertex[1], absDifErr[2], absErr[3], argErr[4], reErr[5], imErr[6], re[7], im[8], gradNormSq[9]");
		ConformalDataList dataList = new ObjectFactory().createConformalDataList();
		ConformalAdapterSet a = new ConformalAdapterSet();
		TypedAdapterSet<double[]> da = a.querySet(double[].class);
		for (int i = 0; i < maxSubdivision; i ++) {
			CoHDS hds = new CoHDS();
			// predefined vertices
			for (int vi = 0; vi < vertices.length; vi++) {
				CoVertex v = hds.addNewVertex();
				v.P = new double[] {vertices[vi][0], vertices[vi][1], vertices[vi][2], 1.0};
				Pn.setToLength(v.P, v.P, 1, Pn.EUCLIDEAN);
			}
			// extra points
			for (int j = 0; j < numExtraPoints; j++) {
				CoVertex v = hds.addNewVertex();
				v.P = new double[] {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian(), 1.0};
				Pn.setToLength(v.P, v.P, 1, Pn.EUCLIDEAN);
			}
			ConvexHull.convexHull(hds, da, 1E-8);
			// subdivision
			for (int si = 0; si < i; si++) {
				CoHDS subdivided = new CoHDS();
				loop.subdivide(hds, subdivided, da);
				hds = subdivided;
				// reset branch data
				for (int vi = 0; vi < vertices.length; vi++) {
					CoVertex v = hds.getVertex(vi);
					v.P = new double[] {vertices[vi][0], vertices[vi][1], vertices[vi][2], 1.0};
					Pn.setToLength(v.P, v.P, 1, Pn.EUCLIDEAN);
				}
				// project to the sphere in every step
				for (CoVertex v : hds.getVertices()) {
					Pn.setToLength(v.P, v.P, 1, Pn.EUCLIDEAN);
				}
			}
			HalfedgeEmbedding input = DataUtility.toHalfedgeEmbedding("input " + i, hds, a, Position4d.class, null);
			dataList.getData().add(input);
			int numVerts = hds.numVertices();
			Complex tau = null;
			CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				DiscreteEllipticUtility.generateEllipticImage(hds, 0, glueSet, branchIndices);
				tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(hds, 1E-9, cutInfo);
			} catch (Exception e) {
				System.out.println("Error: " + e.getMessage());
				continue;
			}
			double absErr = tau.abs() - getExpectedTau().abs();
			double absDifErr = tau.minus(getExpectedTau()).abs();
			double argErr = tau.arg() - getExpectedTau().arg();
			double reErr = tau.re - getExpectedTau().re;
			double imErr = tau.im - getExpectedTau().im;
			writeData(numVerts + "\t" + absDifErr + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + tau.re + "\t" + tau.im + "\t" + EuclideanUnwrapperPETSc.lastGNorm);
			HalfedgeMap map = DataUtility.toHalfedgeMap("map " + i, hds, a, TexturePosition4d.class, Position4d.class, cutInfo);
			dataList.getData().add(map);
		}
		try {
			FileOutputStream fout = new FileOutputStream(fileName + ".xml");
			DataIO.writeConformalDataList(dataList, fout);
		} catch (Exception e) {
			log.log(Level.WARNING, e.getMessage(), e);
		}
	}

}
