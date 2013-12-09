package de.varylab.discreteconformal.plugin.schottky;

import java.io.InputStream;
import java.io.OutputStream;
import java.util.List;

import javax.xml.bind.JAXBException;

import de.varylab.conformallab.data.DataFactory;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.SchottkyData;

public class SchottkyIO {

	
	public static void writeSchottkyData(List<SchottkyGenerator> generators, OutputStream out) throws JAXBException {
		SchottkyData schottkyData = DataUtility.toSchottkyData("Schottky Data", generators);
		DataFactory.writeSchottkyData(schottkyData, out);
	}


	public static List<SchottkyGenerator> readSchottkyData(InputStream in) throws JAXBException {
		SchottkyData data = DataFactory.loadSchottkyData(in);
		List<SchottkyGenerator> result = DataUtility.toGeneratorsList(data);
		return result;
	}

	
}
