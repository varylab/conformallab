package de.varylab.discreteconformal;

import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.Set;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.jreality.util.NativePathUtility;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.calculator.EdgeLengthCalculator;
import de.varylab.discreteconformal.plugin.EllipticModulusEngine;
import de.varylab.discreteconformal.util.Delaunay;

public class Convergence {

	public static void main(String[] args) throws Exception {
		NativePathUtility.set("native");
		
		OptionParser p = new OptionParser();
		OptionSpec<String> fileBase = p.accepts("f", "Name of the data series").withRequiredArg().ofType(String.class);
		OptionSpec<Integer> minPoints = p.accepts("min", "Minimum number of extra points").withRequiredArg().ofType(Integer.class);
		OptionSpec<Integer> maxPoints = p.accepts("max", "Maximum number of extra points").withRequiredArg().ofType(Integer.class);
		OptionSpec<Complex> expectedTau = p.accepts("t").withRequiredArg().ofType(Complex.class);
		OptionSet opts = p.parse(args);
		
		File errFile = new File(fileBase.value(opts) + "_Err.dat");
		FileWriter fwErr = new FileWriter(errFile);
		fwErr.write("# index, absErr, argErr, reErr, imErr");
		for (int i = 1; i <= 5000; i++) {
			System.out.println(i + " extra points --------------------");
			Set<CoEdge> glueSet = new HashSet<CoEdge>();
			Set<CoEdge> cutSet = new HashSet<CoEdge>();
			Complex tau = null;
//			double maxCircRad = 0.0;
			
			CoHDS hds = new CoHDS();
			CoVertex v1 = hds.addNewVertex();
			CoVertex v2 = hds.addNewVertex();
			CoVertex v3 = hds.addNewVertex();
			CoVertex v4 = hds.addNewVertex();
			v1.getPosition().set(1, 0, 0);
			v2.getPosition().set(0, 1, 0);
			v3.getPosition().set(-1, 0, 0);
			v4.getPosition().set(0, -1, 0);
//			double a = 2 * sin(PI / 4);
//			v1.getPosition().set(a, 1, 0);
//			v2.getPosition().set(-a, 1, 0);
//			v3.getPosition().set(0, -1, -a);
//			v4.getPosition().set(0, -1, a);
			
			try {
				EllipticModulusEngine.generateEllipticCurve(hds, i, glueSet, cutSet);
				Delaunay.constructDelaunay(hds, new EdgeLengthCalculator());
				tau = EllipticModulusEngine.calculateModulus(hds);
//				maxCircRad = getMaxCircumRadius(hds);
			} catch (Exception e) {
				continue;
			}
			double absErr = Math.abs(tau.abs() - 1);
			double argErr = Math.abs(tau.arg() - 1.23456);
			double reErr = Math.abs(Math.abs(tau.re) - 0);
			double imErr = Math.abs(Math.abs(tau.im) - 1);
			fwErr.write(i + "\t" + absErr + "\t" + argErr + "\t" + reErr + "\t" + imErr + "\n");
			fwErr.flush();
		}
		fwErr.close();
	}

}
