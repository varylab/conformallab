package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.Clausen.Л;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.calculateTetrahedronVolume;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.calculateTetrahedronVolumeWithIdealVertexAtGamma;
import static de.varylab.discreteconformal.functional.HyperIdealUtility.ζ_13;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static org.junit.Assert.*;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

public class HyperIdealUtilityTest {

	@Test@Ignore
	public void testZeta13() throws Exception {
		double ζ_13 = ζ_13(0.1, 0.1, 0.1);
		double ζ_e = 0.1;
		Assert.assertEquals(ζ_e, ζ_13, 1E-15);
	}
	
	@Test
	public void testVolumeEuclidean() throws Exception {
		double βi = acos(1.0 / 3);
		double βj = acos(1.0 / 3);
		double βk = acos(1.0 / 3);
		double V = calculateTetrahedronVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(0.0, V, 1E-7);
	}
	
	@Test
	public void testVolumeRegularIdeal1() throws Exception {
		double βi = PI / 3;
		double βj = PI / 3;
		double βk = PI / 3;
		double Ve = Л(βi) + Л(βj) + Л(βk);
		double V = calculateTetrahedronVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testVolumeRegularIdeal2() throws Exception {
		double βi = PI / 2;
		double βj = PI / 4;
		double βk = PI / 4;
		double Ve = Л(βi) + Л(βj) + Л(βk);
		double V = calculateTetrahedronVolume(βi, βj, βk, βi, βj, βk);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testVolumeOctahedron() throws Exception {
		double Ve = 8*Л(PI/4);
		double V = calculateTetrahedronVolume(0, 0, 0, 0, 0, 0);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testVolumeSingleHyperidealVertex() throws Exception {
		double βi = PI / 5;
		double βj = PI / 4;
		double βk = PI / 4;
		double ai = (PI + βi - βj - βk) / 2;
		double aj = (PI + βj - βi - βk) / 2;
		double ak = (PI + βk - βi - βj) /2;
		double aijk = (PI - βk - βi - βj) / 2;
		double Ve = 0.5 * (Л(βi) + Л(βj) + Л(βk) + Л(ai) + Л(aj) + Л(ak) + Л(aijk));
		double V = calculateTetrahedronVolume(βi, βj, βk, ai, aj, ak);
		assertEquals(Ve, V, 1E-12);
	}
	
	@Test
	public void testVolumeWithDegenerateTriangle() throws Exception {
		double V = calculateTetrahedronVolume(0.0, PI, 0.0, 0.0, 0.0, PI);
		Assert.assertEquals(0.0, V, 1E-12);
		Assert.assertFalse("Volume must not be NaN", Double.isNaN(V));
	}
	
	@Test
	public void testCompareGeneralAndSingleIdealVolumeFormulas() {
		double EPS = 0;
		double βi = PI / 3;
		double βj = PI / 3;
		double βk = PI / 3;
		double ai = PI / 3 - EPS;
		double aj = PI / 3 - EPS;
		double ak = PI / 3 - EPS;
		double Ve = calculateTetrahedronVolumeWithIdealVertexAtGamma(βi, βj, βk, ai, aj, ak);
		double V = calculateTetrahedronVolume(βi, βj, βk, ai, aj, ak);
		assertEquals(Ve, V, 1E-12);
	}
	
}
