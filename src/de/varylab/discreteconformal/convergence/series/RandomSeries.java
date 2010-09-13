package de.varylab.discreteconformal.convergence.series;

import org.junit.Test;

import de.varylab.discreteconformal.convergence.ConvergenceSeries;

public class RandomSeries {

//	@Test
//	public void pointsMultiRatioExp2() throws Exception {
//		ConvergenceSeries.performConvergenceSeries(
//			"-M", "Random",
//			"--base", "data/convergence/pointsMultiRatioExp2",
//			"--name", "random02-thresh20",
//			"--pin", "data/convergence/random02.obj",
//			"--min", "0",
//			"--max", "2000",
//			"--nopt", "0",
//			"--QM", "MeanMultiRatio",
//			"--QE", "2.0",
//			"--QT", "20.0"
//		);
//	}
	
	@Test
	public void pointsMultiRatioExp10() throws Exception {
		ConvergenceSeries.performConvergenceSeries(
			"-M", "Random",
			"--base", "data/convergence/pointsMultiRatioExp2",
			"--name", "random02",
			"--pin", "data/convergence/random02.obj",
			"--min", "0",
			"--max", "2000",
			"--nopt", "0",
			"--QM", "MeanMultiRatio",
			"--QE", "10.0"
		);
	}
	
}
