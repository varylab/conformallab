package de.varylab.discreteconformal;

import de.jreality.util.NativePathUtility;
import de.varylab.discreteconformal.convergence.ConvergenceSeries;
import de.varylab.discreteconformal.logging.LoggingUtility;

public class Convergence {

	public static void main(String[] args) {
		NativePathUtility.set("native");
		LoggingUtility.initLogging();
		try {
			ConvergenceSeries.performConvergenceSeries(args);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
