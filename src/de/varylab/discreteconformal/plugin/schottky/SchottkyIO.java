package de.varylab.discreteconformal.plugin.schottky;

import java.io.InputStream;
import java.io.OutputStream;
import java.util.LinkedList;
import java.util.List;

import javax.xml.bind.JAXBException;

import de.varylab.conformallab.data.DataFactory;
import de.varylab.conformallab.data.types.Circle;
import de.varylab.conformallab.data.types.ObjectFactory;
import de.varylab.conformallab.data.types.SchottkyData;

public class SchottkyIO {

	
	public static SchottkyData toSchottkyData(String name, List<SchottkyGenerator> generators) {
		ObjectFactory of = new ObjectFactory();
		SchottkyData schottkyData = of.createSchottkyData();
		schottkyData.setName(name);
		for (SchottkyGenerator s : generators) {
			de.varylab.conformallab.data.types.SchottkyGenerator g = of.createSchottkyGenerator();
			de.varylab.conformallab.data.types.Complex Ad = of.createComplex();
			de.varylab.conformallab.data.types.Complex Bd = of.createComplex();
			de.varylab.conformallab.data.types.Complex Mud = of.createComplex();
			de.varylab.conformallab.data.types.Complex center = of.createComplex();
			Circle Cd = of.createCircle();
			Ad.setRe(s.getA().re);
			Ad.setIm(s.getA().im);
			Bd.setRe(s.getB().re);
			Bd.setIm(s.getB().im);
			Mud.setRe(s.getMu().re);
			Mud.setIm(s.getMu().im);
			center.setRe(s.getCycle().getCenter().re);
			center.setIm(s.getCycle().getCenter().im);
			Cd.setCenter(center);
			Cd.setRadius(s.getCycle().getRadius());
			g.setA(Ad);
			g.setB(Bd);
			g.setMu(Mud);
			g.setCircle(Cd);
			schottkyData.getGenerators().add(g);
		}
		return schottkyData;
	}
	
	
	public static void writeSchottkyData(List<SchottkyGenerator> generators, OutputStream out) throws JAXBException {
		SchottkyData schottkyData = toSchottkyData("Schottky Data", generators);
		DataFactory.writeSchottkyData(schottkyData, out);
	}


	public static List<SchottkyGenerator> toGeneratorList(InputStream in) throws JAXBException {
		List<SchottkyGenerator> result = new LinkedList<SchottkyGenerator>();
		SchottkyData data = DataFactory.loadSchottkyData(in);
		for (de.varylab.conformallab.data.types.SchottkyGenerator s : data.getGenerators()) {
			SchottkyGenerator g = new SchottkyGenerator();
			g.getA().setRe(s.getA().getRe());
			g.getA().setIm(s.getA().getIm());
			g.getB().setRe(s.getB().getRe());
			g.getB().setIm(s.getB().getIm());
			g.getMu().setRe(s.getMu().getRe());
			g.getMu().setIm(s.getMu().getIm());
			g.getCycle().getCenter().setRe(s.getCircle().getCenter().getRe());
			g.getCycle().getCenter().setIm(s.getCircle().getCenter().getIm());
			g.getCycle().setRadius(s.getCircle().getRadius());
			result.add(g);
		}
		return result;
	}
	
	
	public static List<SchottkyGenerator> readSchottkyData(InputStream in) throws JAXBException {
		List<SchottkyGenerator> result = toGeneratorList(in);
		return result;
	}

	
}
