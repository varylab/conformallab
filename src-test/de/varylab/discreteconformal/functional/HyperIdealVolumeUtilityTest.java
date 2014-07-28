package de.varylab.discreteconformal.functional;

import static de.varylab.discreteconformal.functional.Clausen.Л;
import static java.lang.Math.PI;

import org.junit.Assert;
import org.junit.Test;

public class HyperIdealVolumeUtilityTest {

	
	@Test
	public void testAgainstLobachewskyVolume() throws Exception {
		double βi = Math.PI / 3 - 0.1;
		double βj = Math.PI / 3 - 0.1;
		double βk = Math.PI / 3 - 0.1;
		double ai = PI + βi - βj - βk;
		double aj = PI + βj - βi - βk;
		double ak = PI + βk - βi - βj;
		double aijk = PI - βk - βi - βj;
		double Vexpected = 0.5 * (Л(βi) + Л(βj) + Л(βk) + Л(ai/2) + Л(aj/2) + Л(ak/2) + Л(aijk/2));
		double V = HyperIdealVolumeUtility.calculateVolume(βi, βj, βk, 0, 0, 0);
		Assert.assertEquals(Vexpected, V, 1E-12);
	}
	
}
