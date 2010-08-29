package de.varylab.discreteconformal;

import static de.jreality.scene.data.Attribute.COORDINATES;
import static java.lang.Math.sqrt;
import geom3d.Point;

import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.PointSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.io.HalfedgeIO;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.plugin.EllipticModulusEngine;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.unwrapper.SphereUtility;

public class Convergence {

	private static Random
		rnd = new Random();
	
	
	public static void main(String[] args) throws Exception {
		NativePathUtility.set("native");
		OptionParser p = new OptionParser();
		OptionSpec<String> fileBaseSpec = p.accepts("f", "Base file name of the data series").withRequiredArg().ofType(String.class);
		OptionSpec<Integer> minPointsSpec = p.accepts("min", "Minimum number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> maxPointsSpec = p.accepts("max", "Maximum number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(50);
		OptionSpec<Integer> incPointsSpec = p.accepts("inc", "Increment").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		OptionSpec<Double> expectedTauReSpec = p.accepts("tre").withRequiredArg().ofType(Double.class).defaultsTo(0.5);
		OptionSpec<Double> expectedTauImSpec = p.accepts("tim").withRequiredArg().ofType(Double.class).defaultsTo(sqrt(3) / 2);
		OptionSpec<String> inputObj = p.accepts("pin", "Predefined branch points as obj file").withRequiredArg().ofType(String.class);
		OptionSpec<Integer> optStepsSpec = p.accepts("nopt", "Number of triangulation optimization steps").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		p.accepts("reuse", "Reuse extra points");
		p.accepts("help", "Prints Help Information");
		OptionSet opts = p.parse(args);
		
		int minExtraPoints = minPointsSpec.value(opts);
		int maxExtraPoints = maxPointsSpec.value(opts);
		int incExtraPoints = incPointsSpec.value(opts);
		
		double tauReExp = opts.valueOf(expectedTauReSpec);
		double tauImExp = opts.valueOf(expectedTauImSpec);
		Complex tauExp = new Complex(tauReExp, tauImExp);
		System.out.println("Expected value if tau: " + tauExp);

		boolean reuse = opts.has("reuse");
		boolean trianopt = opts.has("nopt");
		int numOptSteps = 0;
		if (trianopt) {
			numOptSteps = opts.valueOf(optStepsSpec);
			if (numOptSteps == 0) {
				trianopt = false;
			}
		}
		
		if (opts.has("help")) {
			p.printHelpOn(System.out);
			return;	
		}
		if (!opts.has(fileBaseSpec)) {
			p.printHelpOn(System.out);
			return;
		}
		
		File errFile = new File(fileBaseSpec.value(opts) + "_Err.dat");
		if (errFile.exists()) {
			System.err.println("File " + errFile + " exists. Overwrite? (y/n)");
			char c = (char)System.in.read();
			if (c != 'y') return;
		}
		FileWriter fwErr = new FileWriter(errFile);
		fwErr.write("# index[1], absErr[2], argErr[3], reErr[4], imErr[5], gradNormSq[6]\n");
		for (int i = minExtraPoints; i <= maxExtraPoints; i += incExtraPoints) {
			if (reuse) {
				rnd.setSeed(123);
			}
			System.out.println(i + " extra points --------------------");
			Set<CoEdge> glueSet = new HashSet<CoEdge>();
			Set<CoEdge> cutSet = new HashSet<CoEdge>();
			Complex tau = null;
			CoHDS hds = new CoHDS();
			CoVertex v1 = hds.addNewVertex();
			CoVertex v2 = hds.addNewVertex();
			CoVertex v3 = hds.addNewVertex();
			CoVertex v4 = hds.addNewVertex();
			if (opts.hasArgument(inputObj)) {
				String objIn = opts.valueOf(inputObj);
				ReaderOBJ objReader = new ReaderOBJ();
				SceneGraphComponent c = objReader.read(new File(objIn));
				PointSet g = (PointSet)SceneGraphUtility.getFirstGeometry(c);
				double[][] vertices = g.getVertexAttributes(COORDINATES).toDoubleArrayArray(null);
				if (vertices.length < 4) {
					throw new RuntimeException("Not enough vertices in file " + objIn);
				}
				v1.getPosition().set(vertices[0][0], vertices[0][1], vertices[0][2]);
				v2.getPosition().set(vertices[1][0], vertices[1][1], vertices[1][2]);
				v3.getPosition().set(vertices[2][0], vertices[2][1], vertices[2][2]);
				v4.getPosition().set(vertices[3][0], vertices[3][1], vertices[3][2]);
			} else {
				v1.getPosition().set(1, 1, 1);
				v2.getPosition().set(1, -1, -1);
				v3.getPosition().set(-1, 1, -1);
				v4.getPosition().set(-1, -1, 1);
			}
			// additional points
			for (int j = 0; j < i; j++) {
				double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
				Point pointPos = new Point(pos);
				pointPos.normalize();
				CoVertex v = hds.addNewVertex();
				v.setPosition(pointPos);
			}
			// optimize triangulation
			if (trianopt) {
				Set<CoVertex> fixedVertices = new HashSet<CoVertex>();
				fixedVertices.add(v1);
				fixedVertices.add(v2);
				fixedVertices.add(v3);
				fixedVertices.add(v4);
				SphereUtility.equalizeSphereVertices(hds, fixedVertices, numOptSteps, 1E-6);
				AdapterSet a = new AdapterSet(new PositionAdapter());
				String fileName = "data/convergence/objseries/optGeom" + i + ".obj";
				HalfedgeIO.writeOBJ(hds, a, fileName);
			}
			
			try {
				EllipticModulusEngine.generateEllipticCurve(hds, 0, glueSet, cutSet);
				tau = EllipticModulusEngine.calculateModulus(hds);
			} catch (Exception e) {
				System.err.println(e.getLocalizedMessage());
				continue;
			}
			double absErr = tau.abs() - tauExp.abs();
			double argErr = tau.arg() - tauExp.arg();
			double reErr = tau.re - tauExp.re;
			double imErr = tau.im - tauExp.im;
			System.out.println("tau = " + tau);
			fwErr.write(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + EuclideanUnwrapperPETSc.lastGNorm + "\n");
			fwErr.flush();
		}
		fwErr.close();
	}

}
