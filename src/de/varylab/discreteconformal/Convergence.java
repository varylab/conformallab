package de.varylab.discreteconformal;

import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.Set;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.Geometry;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.PositionAdapter;
import de.varylab.discreteconformal.heds.calculator.EdgeLengthCalculator;
import de.varylab.discreteconformal.plugin.EllipticModulusEngine;
import de.varylab.discreteconformal.util.Delaunay;

public class Convergence {

	public static void main(String[] args) throws Exception {
		NativePathUtility.set("native");
		
		OptionParser p = new OptionParser();
		OptionSpec<String> fileBaseSpec = p.accepts("f", "Name of the data series").withRequiredArg().ofType(String.class);
		OptionSpec<Integer> minPointsSpec = p.accepts("min", "Minimum number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(0);
		OptionSpec<Integer> maxPointsSpec = p.accepts("max", "Maximum number of extra points").withRequiredArg().ofType(Integer.class).defaultsTo(50);
		OptionSpec<Integer> incPointsSpec = p.accepts("inc", "Increment").withRequiredArg().ofType(Integer.class).defaultsTo(1);
//		OptionSpec<Complex> expectedTauSpec = p.accepts("t").withRequiredArg().ofType(Complex.class);
		OptionSpec<String> inputObj = p.accepts("pin", "Predefined branch points as obj file").withRequiredArg().ofType(String.class);
		OptionSet opts = p.parse(args);
		
		int minExtraPoints = minPointsSpec.value(opts);
		int maxExtraPoints = maxPointsSpec.value(opts);
		int incExtraPoints = incPointsSpec.value(opts);
		
//		Complex tauExp = opts.valueOf(expectedTauSpec);
//		System.out.println("Expected value if tau: " + tauExp);
		
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
		fwErr.write("# index, absErr, argErr, reErr, imErr\n");
		for (int i = minExtraPoints; i <= maxExtraPoints; i += incExtraPoints) {
			System.out.println(i + " extra points --------------------");
			Set<CoEdge> glueSet = new HashSet<CoEdge>();
			Set<CoEdge> cutSet = new HashSet<CoEdge>();
			Complex tau = null;
			CoHDS hds = new CoHDS();
			if (opts.hasArgument(inputObj)) {
				String objIn = opts.valueOf(inputObj);
				ReaderOBJ objReader = new ReaderOBJ();
				SceneGraphComponent c = objReader.read(new File(objIn));
				ConverterJR2Heds converter = new ConverterJR2Heds();
				Geometry g = SceneGraphUtility.getFirstGeometry(c);
				if (g instanceof IndexedFaceSet) {
					AdapterSet aSet = new AdapterSet(new PositionAdapter());
					converter.ifs2heds((IndexedFaceSet)g, hds, aSet);
				}
			} else {
				CoVertex v1 = hds.addNewVertex();
				CoVertex v2 = hds.addNewVertex();
				CoVertex v3 = hds.addNewVertex();
				CoVertex v4 = hds.addNewVertex();
				double a = 2 * Math.sin(Math.PI / 4);
				v1.getPosition().set(a, 1, 0);
				v2.getPosition().set(-a, 1, 0);
				v3.getPosition().set(0, -1, -a);
				v4.getPosition().set(0, -1, a);
			}
			try {
				EllipticModulusEngine.generateEllipticCurve(hds, i, glueSet, cutSet);
				Delaunay.constructDelaunay(hds, new EdgeLengthCalculator());
				tau = EllipticModulusEngine.calculateModulus(hds);
			} catch (Exception e) {
				System.err.println(e.getLocalizedMessage());
				continue;
			}
			double absErr = Math.abs(tau.abs() - 1);
			double argErr = Math.abs(tau.arg() - Math.PI/4);
			double reErr = Math.abs(Math.abs(tau.re) - 0.5);
			double imErr = Math.abs(Math.abs(tau.im) - Math.sqrt(3)/3);
			fwErr.write(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\n");
			fwErr.flush();
		}
		fwErr.close();
	}

}
