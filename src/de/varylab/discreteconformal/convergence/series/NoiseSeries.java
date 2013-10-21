package de.varylab.discreteconformal.convergence.series;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;


public class NoiseSeries {

	public void icosahedron011011Reference() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Noise",
			"--base", "data/convergence/noise",
			"--name", "icosahedron-0-1-10-11-Reference",
			"--pin", "data/convergence/icosahedron.obj",
			"--bpi", "0,1,10,11",
			"--IT", "1000",
			"--NC", "0.0"
		);
	}
	
//	@Test
//	public void icosahedron011011() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Noise",
//			"--base", "data/convergence/noise",
//			"--name", "icosahedron-0-1-10-11",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "0,1,10,11",
//			"--IT", "100000",
//			"--NC", "0.05"
//		);
//	}
//	
}
