package de.varylab.discreteconformal;

import static de.jreality.scene.data.Attribute.COORDINATES;
import static java.lang.Math.sqrt;

import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.Set;
import java.util.StringTokenizer;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.PointSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.algorithm.subdivision.LoopLinear;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.calculator.SubdivisionCalculator;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class ConvergenceSubdivision {

	
	/**
	 * Calculates series of elliptic half-period ratios for a given set of
	 * branch points.
	 * It uses loop subdivision to increase the resolution of the given convex hull of
	 * four given branch points.
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		NativePathUtility.set("native");
		OptionParser p = new OptionParser();
		OptionSpec<String> fileBaseSpec = p.accepts("f", "Base file name of the data series").withRequiredArg().ofType(String.class);
		OptionSpec<Integer> mixSubdivisionSpec = p.accepts("max", "Maxmum number of subdivision steps").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Double> expectedTauReSpec = p.accepts("tre").withRequiredArg().ofType(Double.class).defaultsTo(0.5);
		OptionSpec<Double> expectedTauImSpec = p.accepts("tim").withRequiredArg().ofType(Double.class).defaultsTo(sqrt(3) / 2);
		OptionSpec<String> inputObj = p.accepts("pin", "Predefined branch points as obj file").withRequiredArg().ofType(String.class);
		OptionSpec<String> branchIndicesSpec = p.accepts("bpi", "Indices of the four branch points in the obj file").withRequiredArg().ofType(String.class).defaultsTo("0,1,2,3");
		p.accepts("help", "Prints Help Information");
		OptionSet opts = p.parse(args);
		
		int maxSubdivision = mixSubdivisionSpec.value(opts);
		
		double tauReExp = opts.valueOf(expectedTauReSpec);
		double tauImExp = opts.valueOf(expectedTauImSpec);
		Complex tauExp = new Complex(tauReExp, tauImExp);
		System.out.println("Expected value if tau: " + tauExp);

		if (opts.has("help")) {
			p.printHelpOn(System.out);
			return;	
		}
		if (!opts.has(fileBaseSpec) || !opts.has(inputObj)) {
			p.printHelpOn(System.out);
			return;
		}
		
		String fileBase = fileBaseSpec.value(opts);
		File errFile = new File(fileBase + "/err.dat");
		if (errFile.exists()) {
			System.err.println("File " + errFile + " exists. Overwrite? (y/n)");
			char c = (char)System.in.read();
			if (c != 'y') return;
		} else {
			File dirFile = new File(fileBase);
			dirFile.mkdirs();
		}
		
		int[] branchIndices = {0,1,2,3};
		String branchIndexString = branchIndicesSpec.value(opts);
		StringTokenizer bst = new StringTokenizer(branchIndexString, ",");
		int actIndex = 0;
		while (bst.hasMoreTokens()) {
			try {
				int index = Integer.parseInt(bst.nextToken());
				branchIndices[actIndex++] = index;
			} catch (Exception e) {
				e.printStackTrace();
				return;
			}
		}
		
		// get predefined vertices
		String objIn = opts.valueOf(inputObj);
		ReaderOBJ objReader = new ReaderOBJ();
		SceneGraphComponent c = objReader.read(new File(objIn));
		PointSet g = (PointSet)SceneGraphUtility.getFirstGeometry(c);
		double[][] vertices = g.getVertexAttributes(COORDINATES).toDoubleArrayArray(null);
		if (vertices.length < 4) {
			throw new RuntimeException("Not enough vertices in file " + objIn);
		}
		
		LoopLinear loop = new LoopLinear();
		SubdivisionCalculator sc = new SubdivisionCalculator();
		
		FileWriter fwErr = new FileWriter(errFile);
		fwErr.write("# numVertex[1], absErr[2], argErr[3], reErr[4], imErr[5], gradNormSq[6]\n");
		for (int i = 0; i < maxSubdivision; i ++) {
			System.out.println("subdivision step " + i + " --------------------");
			CoHDS hds = new CoHDS();
			// predefined vertices
			for (int vi = 0; vi < vertices.length; vi++) {
				CoVertex v = hds.addNewVertex();
				v.getPosition().set(vertices[vi][0], vertices[vi][1], vertices[vi][2]);	
			}
			ConvexHull.convexHull(hds, sc, 1E-8);
			// subdivision
			for (int si = 0; si < i; si++) {
				CoHDS subdivided = new CoHDS();
				loop.subdivide(hds, subdivided, sc, sc, sc);
				hds = subdivided;
			}
			
			Complex tau = null;
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				DiscreteEllipticUtility.generateEllipticImage(hds, 0, glueSet, branchIndices);
				tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(hds, 1E-8);
			} catch (Exception e) {
				System.out.println("Error: " + e.getMessage());
				continue;
			}
			double absErr = tau.abs() - tauExp.abs();
			double argErr = tau.arg() - tauExp.arg();
			double reErr = tau.re - tauExp.re;
			double imErr = tau.im - tauExp.im;
			System.out.println("tau = " + tau);
			String writeLine = hds.numVertices() + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + EuclideanUnwrapperPETSc.lastGNorm;
			System.out.println(writeLine);
			fwErr.write(writeLine + "\n");
			fwErr.flush();
		}
		fwErr.close();
	}

}
