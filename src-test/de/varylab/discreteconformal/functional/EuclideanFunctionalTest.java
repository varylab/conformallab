package de.varylab.discreteconformal.functional;

import static de.jreality.util.SceneGraphUtility.getFirstGeometry;

import java.io.IOException;
import java.util.Random;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class EuclideanFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	private CTheta
		theta = new CTheta();
	private CVariable
		variable = new CVariable();
	private CLambda
		lambda = new CLambda();
	private CAlpha
		alpha = new CAlpha();
	private CInitialEnergy
		energy = new CInitialEnergy();
	public EuclideanFunctional<CoVertex, CoEdge, CoFace>
		functional = new EuclideanFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	
	@Override
	public void init() {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		CoHDS hds = new CoHDS(); 
		try {
			Input in = new Input("Obj File", EuclideanFunctionalTest.class.getResourceAsStream("square01.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)getFirstGeometry(c);
			ConverterJR2Heds converter = new ConverterJR2Heds();
			hds = new CoHDS();
			converter.ifs2heds(ifs, hds, a, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		int n = UnwrapUtility.prepareInvariantDataEuclidean(functional, hds, a);
		Random rnd = new Random(); 
		rnd.setSeed(1);
		
		Vector x = new DenseVector(n);
		for (Integer i = 0; i < x.size(); i++) {
			x.set(i, rnd.nextDouble() - 0.5);
		}
		MyDomainValue u = new MyDomainValue(x);
		
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
		setXHessian(u);
	}
	
	
}
