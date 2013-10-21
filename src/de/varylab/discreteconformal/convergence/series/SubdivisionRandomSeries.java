package de.varylab.discreteconformal.convergence.series;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class SubdivisionRandomSeries {

	
//	@Test
//	public void random02() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Subdivision",
//			"--base", "data/convergence/subdivision",
//			"--name", "random02",
//			"--pin", "data/convergence/random02.obj",
//			"--bpi", "0,1,2,3",
//			"--max", "6"
//		);
//	}
//	
//	@Test
//	public void random03() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Subdivision",
//			"--base", "data/convergence/subdivision",
//			"--name", "random03",
//			"--pin", "data/convergence/random03.obj",
//			"--bpi", "0,1,2,3",
//			"--max", "6"
//		);
//	}
//	
//	@Test
//	public void random04() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Subdivision",
//			"--base", "data/convergence/subdivision",
//			"--name", "random04",
//			"--pin", "data/convergence/random04.obj",
//			"--bpi", "0,1,2,3",
//			"--max", "6"
//		);
//	}
	
	
	public void random02_extra4() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Subdivision",
			"--base", "data/convergence/subdivision",
			"--name", "random02_extra4",
			"--pin", "data/convergence/random02.obj",
			"--bpi", "0,1,2,3",
			"--max", "6",
			"--extra", "4"
		);
	}
	
	public void random03_extra4() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Subdivision",
			"--base", "data/convergence/subdivision",
			"--name", "random03_extra4",
			"--pin", "data/convergence/random03.obj",
			"--bpi", "0,1,2,3",
			"--max", "6",
			"--extra", "4"
		);
	}
	
	public void random04_extra4() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Subdivision",
			"--base", "data/convergence/subdivision",
			"--name", "random04_extra4",
			"--pin", "data/convergence/random04.obj",
			"--bpi", "0,1,2,3",
			"--max", "6",
			"--extra", "4"
		);
	}
	
}
