package de.varylab.discreteconformal.data;

import java.io.InputStream;

import org.junit.Test;

import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.UniformizationData;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;

public class UniformizationDataTest {

	@Test
	public void testLoadUniformizationData() throws Exception {
		InputStream in = getClass().getResourceAsStream("LetterB01_canonical_group.xml");
		UniformizationData data = DataIO.readConformalData(UniformizationData.class, in);
		FundamentalPolygon P = DataUtility.toFundamentalPolygon(data);
		System.out.println(P.getDiscreteGroup());
	}
	
}
