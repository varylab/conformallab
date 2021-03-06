package de.varylab.discreteconformal.plugin.schottky;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.LineNumberReader;
import java.io.OutputStream;
import java.util.LinkedList;
import java.util.List;

import org.custommonkey.xmlunit.XMLAssert;
import org.junit.Assert;
import org.junit.Test;

import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.SchottkyData;

public class SchottkyIOTest {

	@Test
	public void testReadSchottkyData() throws Exception {
		InputStream in = getClass().getResourceAsStream("schottkydata.xml");
		SchottkyData sd = DataIO.readConformalData(SchottkyData.class, in);
		List<SchottkyGenerator> data = DataUtility.toSchottkyGeneratorsList(sd);
		Assert.assertEquals(1, data.size());
		SchottkyGenerator s = data.get(0);
		Assert.assertEquals(-1.0, s.getA().re, 1E-24);
		Assert.assertEquals(2.0, s.getA().im, 1E-24);
		Assert.assertEquals(1.0, s.getB().re, 1E-24);
		Assert.assertEquals(3.0, s.getB().im, 1E-24);
		Assert.assertEquals(0.5, s.getCycle().getRadius(), 1E-24);
		Assert.assertEquals(0.5, s.getCycle().getRadius(), 1E-24);
		Assert.assertEquals(-1.0, s.getCycle().getCenter().re, 1E-24);
		Assert.assertEquals(0.0, s.getCycle().getCenter().im, 1E-24);
	}
	
	@Test
	public void testWriteSchottkyData() throws Exception {
		List<SchottkyGenerator> data = new LinkedList<SchottkyGenerator>();
		data.add(new SchottkyGenerator());
		File testFile = File.createTempFile("testWriteSchottkyData", ".xml");
		testFile.deleteOnExit();
		OutputStream out = new FileOutputStream(testFile);
		SchottkyData sd = DataUtility.toSchottkyData("Schottky Data", data);
		DataIO.writeConformalData(sd, out);
		LineNumberReader r = new LineNumberReader(new FileReader(testFile));
		StringBuffer xml = new StringBuffer();
		String line = r.readLine();
		while (line != null) {
			xml.append(line + "\n");
			line = r.readLine();
		}
		r.close();
		XMLAssert.assertXMLEqual(
			"<SchottkyData name='Schottky Data' xmlns='http://www.varylab.com/conformallab/types'>\n" + 
			"    <SchottkyGenerator>\n" + 
			"        <A im='0.0' re='-1.0'/>\n" + 
			"        <B im='0.0' re='1.0'/>\n" + 
			"        <Mu im='0.0' re='0.05'/>\n" + 
			"        <Circle radius='0.5'>\n" + 
			"            <Center im='0.0' re='-1.0'/>\n" + 
			"        </Circle>\n" + 
			"    </SchottkyGenerator>\n" + 
			"</SchottkyData>", 
			xml.toString()
		);
	}


}
