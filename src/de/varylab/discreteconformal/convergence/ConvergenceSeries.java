package de.varylab.discreteconformal.convergence;

import static de.jreality.scene.data.Attribute.COORDINATES;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;

import joptsimple.OptionParser;
import joptsimple.OptionSet;
import joptsimple.OptionSpec;

import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkFactory;

import de.jreality.math.Pn;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.PointSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.NativePathUtility;
import de.jreality.util.SceneGraphUtility;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public abstract class ConvergenceSeries {

	
	public static enum SeriesMethod {
		Quality,
		Random,
		Subdivision,
		Noise
	}
	
	private static Logger 
		log = Logger.getLogger(ConvergenceSeries.class.getName());
	protected Random
		rnd = new Random();
	protected double[][]
	    vertices = {};
	protected int[]
	    branchIndices = {0,1,2,3};
	protected FileWriter
		errorWriter = null;
	protected Complex
		tauExpected = new Complex();
	
	
	public ConvergenceSeries() {
	}
	
	public ConvergenceSeries(
		double[][] vertices, 
		int[] branchIndices,
		Complex tauExpected,
		FileWriter errorWriter
	) {
		super();
		this.vertices = vertices;
		this.branchIndices = branchIndices;
		this.errorWriter = errorWriter;
		this.tauExpected = tauExpected;
	}


	public static void performConvergenceSeries(String... args) throws Exception {
		NativePathUtility.set("native");
		boolean mFound = false;
		SeriesMethod method = null;
		for (String arg : args) {
			if (mFound) {
				method = SeriesMethod.valueOf(arg);
				break;
			}
			if (arg.equals("-M")) mFound = true;
		}
		if (method == null) {
			throw new RuntimeException("Unknown series method");
		}
		
		ConvergenceSeries series = new ConvergenceQuality();
		switch (method) {
		case Quality:
			series = new ConvergenceQuality();
			break;
		case Random:
			series = new ConvergenceRandom();
			break;
		case Subdivision:
			series = new ConvergenceSubdivision();
			break;
		case Noise:
			series = new ConvergenceNoise();
			break;
		}
		assert series != null;
		
		OptionParser p = new OptionParser();
		OptionSpec<String> methodSpec = p.accepts("M", "Convergence series method: Quality | Random | Subdivision").withRequiredArg().ofType(String.class);
		OptionSpec<String> fileBaseSpec = p.accepts("base", "Base directory of the data series").withRequiredArg().ofType(String.class);
		OptionSpec<String> fileNameSpec = p.accepts("name", "Name of the data series").withRequiredArg().ofType(String.class);
		OptionSpec<String> inputObj = p.accepts("pin", "Predefined branch points as obj file").withRequiredArg().ofType(String.class);
		OptionSpec<String> branchIndicesSpec = p.accepts("bpi", "Indices of the four branch points in the obj file").withRequiredArg().ofType(String.class).defaultsTo("0,1,2,3");
		OptionSpec<Double> tauReSpec = p.accepts("tre", "Real part of expected \tau").withOptionalArg().ofType(Double.class);
		OptionSpec<Double> tauImSpec = p.accepts("tim", "Imaginary part of expected \tau").withOptionalArg().ofType(Double.class);
		p.accepts("help", "Prints Help Information");
		OptionSet opts = series.configureAndParseOptions(p, args);
		
		// print help if options error
		if (opts.has("help") || !opts.has(fileBaseSpec) || !opts.has(fileNameSpec) || !opts.has(methodSpec)) {
			p.printHelpOn(System.out);
			return;
		}
		
		// get file base directory
		String fileBase = fileBaseSpec.value(opts);
		String fileName = fileNameSpec.value(opts);
		File errFile = new File(fileBase + "/" + fileName + ".dat");
		if (errFile.exists()) {
			System.err.println("File " + errFile + " exists. Overwrite? (y/n)");
			char c = (char)System.in.read();
			if (c != 'y') {
				throw new RuntimeException("Exit.");
			}
		} else {
			File dirFile = new File(fileBase);
			dirFile.mkdirs();
		}
		series.errorWriter = new FileWriter(errFile);
		series.writeComment("Series class: " + series.getClass());
		series.writeComment("Error output file " + errFile);
		
		// read vertices
		String objIn = opts.valueOf(inputObj);
		ReaderOBJ objReader = new ReaderOBJ();
		SceneGraphComponent c = objReader.read(new File(objIn));
		PointSet g = (PointSet)SceneGraphUtility.getFirstGeometry(c);
		series.vertices = g.getVertexAttributes(COORDINATES).toDoubleArrayArray(null);
		if (series.vertices.length < 4) {
			throw new RuntimeException("Not enough vertices in file " + objIn);
		}
		series.writeComment("Vertex data from " + objIn);
		
		// branch indices
		series.branchIndices = new int[]{0,1,2,3};
		String branchIndexString = branchIndicesSpec.value(opts);
		StringTokenizer bst = new StringTokenizer(branchIndexString, ",");
		int actIndex = 0;
		while (bst.hasMoreTokens()) {
			try {
				int index = Integer.parseInt(bst.nextToken());
				series.branchIndices[actIndex++] = index;
			} catch (Exception e) {
				throw new RuntimeException("Error while parsing branch point indices: \n\t" + e.getMessage());
			}
		}
		series.writeComment("Branch indices: " + Arrays.toString(series.branchIndices));
		
		
		// calculate half-period ratio with Mathematica
		String[] mlargs = new String[] {
			"-linkmode", "launch", 
			"-linkname", "\"C:\\Program Files\\Wolfram Research\\Mathematica\\7.0\\MathKernel.exe\" " + 
			"-mathlink"
		};
		KernelLink link = null;
		try {
			link = MathLinkFactory.createKernelLink(mlargs);
			link.discardAnswer();
		} catch (Throwable e) {
			log.log(Level.WARNING, "could not create mathematica link: ", e);
		}
		double[] p1 = series.vertices[series.branchIndices[0]];
		double[] p2 = series.vertices[series.branchIndices[1]];
		double[] p3 = series.vertices[series.branchIndices[2]];
		double[] p4 = series.vertices[series.branchIndices[3]];
		if (p1.length == 3) {
			p1 = Pn.homogenize(null, p1);
			p2 = Pn.homogenize(null, p2);
			p3 = Pn.homogenize(null, p3);
			p4 = Pn.homogenize(null, p4);
		}
		Double tauRe = tauReSpec.value(opts);
		Double tauIm = tauImSpec.value(opts);
		if (link != null) {
			series.tauExpected = DiscreteEllipticUtility.calculateHalfPeriodRatioMathLink(p1, p2, p3, p4, link);
			link.close();
		} else if (tauRe != null || tauIm != null) {
			series.tauExpected = new Complex(tauRe, tauIm);
		} else {
			throw new Exception("Could not determine expected tau");
		}
		series.writeComment("Expected tau: " + series.tauExpected);
		
		series.perform();
		series.getErrorWriter().close();
	}
	
	public Random getRnd() {
		return rnd;
	}
	public double[][] getVertices() {
		return vertices;
	}
	public int[] getBranchIndices() {
		return branchIndices;
	}
	public FileWriter getErrorWriter() {
		return errorWriter;
	}
	public Complex getExpectedTau() {
		return tauExpected;
	}
	
	
	protected abstract void perform() throws Exception;
	
	protected OptionSet configureAndParseOptions(OptionParser p, String... args) {
		return p.parse(args);
	}
	
	protected void writeComment(String comment) throws IOException {
		writeErrorLine("#### " + comment + " ####");
	}
	
	
	protected void writeErrorLine(String line) throws IOException {
		log.info(line);
		errorWriter.write(line + "\n");
		errorWriter.flush();
	}
	
}
