package de.varylab.discreteconformal;

import static de.jreality.scene.data.Attribute.COORDINATES;
import static java.lang.Math.max;
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
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.io.HalfedgeIO;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.calculator.SubdivisionCalculator;
import de.varylab.discreteconformal.unwrapper.SphereUtility;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class ConvergenceQuality {

	private static Random
		rnd = new Random();
	
	public static void main(String[] args) throws Exception {
		NativePathUtility.set("native");
		OptionParser p = new OptionParser();
		OptionSpec<String> fileBaseSpec = p.accepts("f", "Base file name of the data series").withRequiredArg().ofType(String.class);
		OptionSpec<Integer> numSamplesSpec = p.accepts("num", "Number of samples").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		OptionSpec<Integer> numExtraPointsSpec = p.accepts("extra", "Number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Double> expectedTauReSpec = p.accepts("tre").withRequiredArg().ofType(Double.class).defaultsTo(0.5);
		OptionSpec<Double> expectedTauImSpec = p.accepts("tim").withRequiredArg().ofType(Double.class).defaultsTo(sqrt(3) / 2);
		OptionSpec<String> inputObj = p.accepts("pin", "Predefined branch points as obj file").withRequiredArg().ofType(String.class);
		OptionSpec<Integer> optStepsSpec = p.accepts("nopt", "Number of triangulation optimization steps").withRequiredArg().ofType(Integer.class).defaultsTo(1);
		p.accepts("reuse", "Reuse extra points");
		p.accepts("help", "Prints Help Information");
		p.accepts("w", "Write OBJ files to objseries folder");
		OptionSet opts = p.parse(args);
		
		int numSamples = numSamplesSpec.value(opts);
		int numExtraPoints = numExtraPointsSpec.value(opts);
		
		double tauReExp = opts.valueOf(expectedTauReSpec);
		double tauImExp = opts.valueOf(expectedTauImSpec);
		Complex tauExp = new Complex(tauReExp, tauImExp);
		System.out.println("Expected value if tau: " + tauExp);

		boolean reuse = opts.has("reuse");
		boolean trianopt = opts.has("nopt");
		boolean writeOBJs = opts.has("w");
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
		
		FileWriter fwErr = new FileWriter(errFile);
		fwErr.write("# index[1], absErr[2], argErr[3], reErr[4], imErr[5], qFun[6]\n");
		for (int i = 0; i <numSamples; i ++) {
			if (reuse) {
				rnd.setSeed(123);
			}
			System.out.println("Taking sample number " + i);
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
			// extra points
			for (int j = 0; j < numExtraPoints; j++) {
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
				if (writeOBJs) {
					File dirFile = new File(fileBase + "/objseries");
					dirFile.mkdirs();
					String fileName = fileBase + "/objseries/optGeom" + i + ".obj";
					try {
						ConvexHull.convexHull(hds, new SubdivisionCalculator(), 1E-8);
					} catch (Exception e) {
						System.err.println("No convex hull for " + fileName);
					}
					HalfedgeIO.writeOBJ(hds, a, fileName);
				}
			}
			try {
				Set<CoEdge> glueSet = new HashSet<CoEdge>();
				DiscreteEllipticUtility.generateEllipticImage(hds, 0, glueSet);
				tau = DiscreteEllipticUtility.calculateHalfPeriodRatio(hds, 1E-8);
			} catch (Exception e) {
				System.err.println(e.getLocalizedMessage());
				continue;
			}
			double absErr = tau.abs() - tauExp.abs();
			double argErr = tau.arg() - tauExp.arg();
			double reErr = tau.re - tauExp.re;
			double imErr = tau.im - tauExp.im;
			double qFun = qFun(hds);
			System.out.println("tau = " + tau);
			System.out.println("qfun = " + qFun);
			fwErr.write(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\t" + qFun + "\n");
			fwErr.flush();
		}
		fwErr.close();
	}
	
	
	
	private static double qFun(CoHDS hds) {
		double r = 0.0;
		for (CoEdge e : hds.getPositiveEdges()) {
			double a = e.getNextEdge().getLength();
			double b = e.getPreviousEdge().getLength();
			double c = e.getOppositeEdge().getNextEdge().getLength();
			double d = e.getOppositeEdge().getPreviousEdge().getLength();
			double q = (a * c) / (b * d);
			double qfun = (q + 1/q)/2 - 1;
			r = max(r, qfun);
		}
		return r;
	}
	
	
	

}