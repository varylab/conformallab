package de.varylab.discreteconformal.util;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.StringReader;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import de.jreality.geometry.IndexedFaceSetUtility;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.functional.DomainValue;
import de.jtem.halfedgetools.functional.Functional;
import de.jtem.halfedgetools.functional.Gradient;
import de.jtem.halfedgetools.functional.Hessian;
import de.jtem.halfedgetools.functional.MyEnergy;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.jtem.jpetsc.Mat;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.TaoAppAddCombinedObjectiveAndGrad;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;

public class TestUtility {

	private static final double
		eps = 1E-5;
	
	public static double[][] parseTaoTestOutput(String output, int len) {
		int fdIndex = output.indexOf("Finite difference gradient");
		int hcIndex = output.indexOf("Hand-coded gradient");
		if (fdIndex == -1 || hcIndex == -1) {
			throw new IllegalArgumentException("Invalid Tao test result.");
		}
		String fdOut = output.substring(fdIndex);
		String hcOut = output.substring(hcIndex);
		LineNumberReader fdReader = new LineNumberReader(new StringReader(fdOut));
		LineNumberReader hcReader = new LineNumberReader(new StringReader(hcOut));
		double[][] result = new double[2][len];
		try {
			String fdNumber = fdReader.readLine();
			String hcNumber = hcReader.readLine();
			int index = 0;
			do {
				fdNumber = fdReader.readLine();
				hcNumber = hcReader.readLine();
				result[0][index] = Double.parseDouble(fdNumber);
				result[1][index++] = Double.parseDouble(hcNumber);
			} while (fdNumber != null && hcNumber != null);
		} catch (Exception e) {}
		return result;
	}
	
	
	public static void calculateFDGradient(TaoAppAddCombinedObjectiveAndGrad app, Vec x, Vec G) {
		for (int i = 0; i < x.getSize(); i++){
			double xi = x.getValue(i);
			x.setValue(i, xi + eps, INSERT_VALUES);
			double f1 = app.evaluateObjectiveAndGradient(x, null);
			x.setValue(i, xi - eps, INSERT_VALUES);
			double f2 = app.evaluateObjectiveAndGradient(x, null);
			double fdGrad = (f1 - f2) / (2 * eps);
			G.setValue(i, fdGrad, INSERT_VALUES);
			x.setValue(i, xi, INSERT_VALUES);
		}
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void calculateFDGradient(HDS hds, Functional<V, E, F> app, int dim, DomainValue x, Gradient G) {
		MyEnergy E = new MyEnergy();
		for (int i = 0; i < dim; i++){
			double xi = x.get(i);
			x.set(i, xi + eps);
			app.evaluate(hds, x, E, null, null);
			double f1 = E.get();
			x.set(i, xi - eps);
			app.evaluate(hds, x, E, null, null);
			double f2 = E.get();
			double fdGrad = (f1 - f2) / (2 * eps);
			G.set(i, fdGrad);
			x.set(i, xi);
		}
	}
	
	public static void calculateFDHessian(TaoAppAddCombinedObjectiveAndGrad app, Vec x, Mat H) {
		double y = app.evaluateObjectiveAndGradient(x, null);
		for (int i = 0; i < x.getSize(); i++){
			for (int j = 0; j < x.getSize(); j++){
				double fdHessian = 0.0;
				double xi = x.getValue(i);
				double xj = x.getValue(j);
				if (i == j) {
					x.setValue(i, xi + eps, INSERT_VALUES);
					double iPlus = app.evaluateObjectiveAndGradient(x, null);
					x.setValue(i, xi + 2*eps, INSERT_VALUES);
					double i2Plus = app.evaluateObjectiveAndGradient(x, null);	
					x.setValue(i, xi - eps, INSERT_VALUES);
					double iMinus = app.evaluateObjectiveAndGradient(x, null);
					x.setValue(i, xi - 2*eps, INSERT_VALUES);
					double i2Minus = app.evaluateObjectiveAndGradient(x, null);
					fdHessian = (-i2Plus/eps + 16*iPlus/eps - 30*y/eps + 16*iMinus/eps - i2Minus/eps) / (12 * eps);
				} else {
					x.setValue(i, xi + eps, INSERT_VALUES);
					x.setValue(j, xj + eps, INSERT_VALUES);
					double iPlusjPlus = app.evaluateObjectiveAndGradient(x, null);
					x.setValue(i, xi + eps, INSERT_VALUES);
					x.setValue(j, xj - eps, INSERT_VALUES);
					double iPlusjMinus = app.evaluateObjectiveAndGradient(x, null);
					x.setValue(i, xi - eps, INSERT_VALUES);
					x.setValue(j, xj + eps, INSERT_VALUES);
					double iMinusjPlus = app.evaluateObjectiveAndGradient(x, null);
					x.setValue(i, xi - eps, INSERT_VALUES);
					x.setValue(j, xj - eps, INSERT_VALUES);
					double iMinusjMinus = app.evaluateObjectiveAndGradient(x, null);
					fdHessian = (iPlusjPlus/eps - iPlusjMinus/eps - iMinusjPlus/eps + iMinusjMinus/eps) / (4 * eps);
				}
				x.setValue(i, xi, INSERT_VALUES);
				x.setValue(j, xj, INSERT_VALUES);
				H.setValue(i, j, fdHessian, INSERT_VALUES);
			}
		}
		H.assemble();
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void calculateFDHessian(HDS hds, Functional<V, E, F> app, int dim, DomainValue x, Hessian H) {
		MyEnergy E = new MyEnergy();
		app.evaluate(hds, x, E, null, null);
		double y = E.get();
		for (int i = 0; i < dim; i++){
			for (int j = 0; j < dim; j++){
				double fdHessian = 0.0;
				double xi = x.get(i);
				double xj = x.get(j);
				if (i == j) {
					x.set(i, xi + eps);
					app.evaluate(hds, x, E, null, null);
					double iPlus = E.get();
					x.set(i, xi + 2*eps);
					app.evaluate(hds, x, E, null, null);
					double i2Plus = E.get();
					x.set(i, xi - eps);
					app.evaluate(hds, x, E, null, null);
					double iMinus = E.get();
					x.set(i, xi - 2*eps);
					app.evaluate(hds, x, E, null, null);
					double i2Minus = E.get();
					fdHessian = (-i2Plus/eps + 16*iPlus/eps - 30*y/eps + 16*iMinus/eps - i2Minus/eps) / (12 * eps);
				} else {
					x.set(i, xi + eps);
					x.set(j, xj + eps);
					app.evaluate(hds, x, E, null, null);
					double iPlusjPlus = E.get();
					x.set(i, xi + eps);
					x.set(j, xj - eps);
					app.evaluate(hds, x, E, null, null);
					double iPlusjMinus = E.get();
					x.set(i, xi - eps);
					x.set(j, xj + eps);
					app.evaluate(hds, x, E, null, null);
					double iMinusjPlus = E.get();
					x.set(i, xi - eps);
					x.set(j, xj - eps);
					app.evaluate(hds, x, E, null, null);
					double iMinusjMinus = E.get();
					fdHessian = (iPlusjPlus/eps - iPlusjMinus/eps - iMinusjPlus/eps + iMinusjMinus/eps) / (4 * eps);
				}
				x.set(i, xi);
				x.set(j, xj);
				H.set(i, j, fdHessian);
			}
		}
	}
	
	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void completeOpposites(Map<E, Double> alphaMap) {
		Set<E> edgeKeys = new HashSet<E>(alphaMap.keySet());
		for (E e : edgeKeys) {
			if (!alphaMap.containsKey(e.getOppositeEdge())) {
				Double val = alphaMap.get(e);
				alphaMap.put(e.getOppositeEdge(), val);
			}
		}
	}
	
	
	public static CoHDS readOBJ(Class<?> root, String name) {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		try {
			Input in = new Input("Obj File", root.getResourceAsStream(name));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			IndexedFaceSetUtility.makeConsistentOrientation(ifs);
			ConverterJR2Heds converter = new ConverterJR2Heds();
			CoHDS hds = new CoHDS();
			AdapterSet a = new AdapterSet();
			a.add(new CoPositionAdapter());
			a.add(new CoTexturePositionAdapter());
			converter.ifs2heds(ifs, hds, a, null);
			return hds;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	
	
	
}
