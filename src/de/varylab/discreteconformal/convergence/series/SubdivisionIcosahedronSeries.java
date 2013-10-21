package de.varylab.discreteconformal.convergence.series;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class SubdivisionIcosahedronSeries {
//
//	@Test
//	public void icosahedron14711() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Subdivision",
//			"--base", "data/convergence/subdivision",
//			"--name", "icosahedron-1-4-7-11",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "1,4,7,11",
//			"--max", "6"
//		);
//	}
//	
//	@Test
//	public void icosahedron15711() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Subdivision",
//			"--base", "data/convergence/subdivision",
//			"--name", "icosahedron-1-5-7-11",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "1,5,7,11",
//			"--max", "6"
//		);
//	}
//	
//	@Test
//	public void icosahedron13711() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Subdivision",
//			"--base", "data/convergence/subdivision",
//			"--name", "icosahedron-1-3-7-11",
//			"--pin", "data/convergence/icosahedron.obj",
//			"--bpi", "1,3,7,11",
//			"--max", "6"
//		);
//	}
	

	public void icosahedron14711Rand() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Subdivision",
			"--base", "data/convergence/subdivision",
			"--name", "icosahedron-1-4-7-11Rand",
			"--pin", "data/convergence/icosahedronRand.obj",
			"--bpi", "1,4,7,11",
			"--max", "6"
		);
	}
	
	public void icosahedron15711() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Subdivision",
			"--base", "data/convergence/subdivision",
			"--name", "icosahedron-1-5-7-11Rand",
			"--pin", "data/convergence/icosahedronRand.obj",
			"--bpi", "1,5,7,11",
			"--max", "6"
		);
	}
	
	public void icosahedron13711() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Subdivision",
			"--base", "data/convergence/subdivision",
			"--name", "icosahedron-1-3-7-11Rand",
			"--pin", "data/convergence/icosahedronRand.obj",
			"--bpi", "1,3,7,11",
			"--max", "6"
		);
	}
	
}
