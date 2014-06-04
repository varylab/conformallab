package de.varylab.discreteconformal;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.logging.Logger;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMap;
import de.varylab.discreteconformal.plugin.TargetGeometry;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class ConformalLabBatch {

	private Logger
		log = Logger.getLogger(ConformalLabBatch.class.getName());
	
	private OptionSpec<TargetGeometry> 
		geometryArg = null;
	private OptionSpec<File>
		fileArg = null;
	
	public void process(String... args) throws Exception {
		OptionParser p = new OptionParser();
		p.accepts("help", "Prints Help Information");
		fileArg = p
			.accepts("file")
			.withRequiredArg()
			.ofType(File.class);
		geometryArg = p
			.accepts("geometry", "Target Geometry")
			.withRequiredArg()
			.ofType(TargetGeometry.class)
			.defaultsTo(TargetGeometry.Automatic);
		OptionSet opts = null;
		try {
			opts = p.parse(args);
			process(opts);
		} catch (Exception e) {
			printUsage(p);
			return;
		}
	}
	
	protected void process(OptionSet opts) throws Exception {
		TargetGeometry tg = geometryArg.value(opts);
		log.info("using target geometry: " + tg);
		InputStream dataIn = null;
		if (opts.hasArgument(fileArg)) {
			dataIn = new FileInputStream(fileArg.value(opts));
		} else {
			dataIn = System.in;
		}
		
		DiscreteMap mapData = DataIO.readConformalData(DiscreteMap.class, dataIn);
		log.info("using conformal data: " + mapData.getName());
		
		DiscreteEmbedding image = mapData.getImage();
//		DiscreteEmbedding domain = mapData.getDomain();
		int genus = DataUtility.calculateGenus(image);
		log.info("loading discreet map of a genus " + genus + " surface");
//		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
//		CoHDS hds = DataUtility.toHDS(image, cutInfo);
		
		if (tg == TargetGeometry.Automatic) {
			tg = UnwrapUtility.calculateTargetGeometry(genus, 0);
			log.info("automatic target geometry: " + tg);
		}
		dataIn.close();
	}
	
	
	private void printUsage(OptionParser p) throws Exception {
		System.out.println("Usage: java -jar DiscreteConformalLab.jar [Options] <data>");
		p.printHelpOn(System.out);
	}
	
}
