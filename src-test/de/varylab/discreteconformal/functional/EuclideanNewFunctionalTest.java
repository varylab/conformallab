package de.varylab.discreteconformal.functional;

import java.io.IOException;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.functional.FunctionalTest;
import de.jtem.halfedgetools.functional.MyDomainValue;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.functional.hds.MyConformalAdapters.CPhi;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomEdgeInfo;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class EuclideanNewFunctionalTest extends FunctionalTest<CoVertex, CoEdge, CoFace> {

	private CTheta
		theta = new CTheta();
	private CPhi
		phi = new CPhi();	
	private CVariable
		variable = new CVariable();
	private CLambda
		lambda = new CLambda();
	private CAlpha
		alpha = new CAlpha();
	private CInitialEnergy
		energy = new CInitialEnergy();
	public EuclideanNewFunctional<CoVertex, CoEdge, CoFace>
		functional = new EuclideanNewFunctional<CoVertex, CoEdge, CoFace>(variable, theta, phi, lambda, alpha, energy);
	
	@Override
	public void init() {
		ReaderOBJ reader = new ReaderOBJ();
		SceneGraphComponent c = null;
		IndexedFaceSet ifs = null;
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		CoHDS hds = new CoHDS(); 
		try {
			Input in = new Input("Obj File", EuclideanNewFunctionalTest.class.getResourceAsStream("square01.obj"));
			c =reader.read(in);
			ifs = (IndexedFaceSet)c.getChildComponent(0).getGeometry();
			ConverterJR2Heds converter = new ConverterJR2Heds();
			hds = new CoHDS();
			converter.ifs2heds(ifs, hds, a, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// one triangle of edges is circular
		for (CoFace f : hds.getFaces()) {
			if (!HalfEdgeUtils.isInteriorFace(f)) continue;
			CoEdge e1 = f.getBoundaryEdge();
			CoEdge e2 = e1.getNextEdge();
			CoEdge e3 = e2.getNextEdge();
			CustomEdgeInfo info = new CustomEdgeInfo();
			info.circularHoleEdge = true;
			e1.info = new CustomEdgeInfo(info);
			e2.info = new CustomEdgeInfo(info);
			e3.info = new CustomEdgeInfo(info);
			e1.getOppositeEdge().info = e1.info;
			e2.getOppositeEdge().info = e2.info;
			e3.getOppositeEdge().info = e3.info;
//			e1.info.phi = Math.PI - 0.3; // one with a modified opposite angle sum
			break;
		}
		
		int n = UnwrapUtility.prepareInvariantDataEuclidean(functional, hds, a);
		Vector x = new DenseVector(n);
		for (CoEdge e : hds.getPositiveEdges()) {
			if (e.getSolverIndex() >= 0) {
				x.set(e.getSolverIndex(), e.getLambda());
			}
		}
		MyDomainValue u = new MyDomainValue(x);
		
		setFunctional(functional);
		setHDS(hds);
		setXGradient(u);
		setXHessian(u);
	}
	
}
