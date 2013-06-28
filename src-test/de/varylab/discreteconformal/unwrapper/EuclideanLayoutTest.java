package de.varylab.discreteconformal.unwrapper;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.SparseVector;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.reader.ReaderOBJ;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.util.Input;
import de.jreality.util.SceneGraphUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.io.HalfedgeIO;
import de.jtem.halfedgetools.jreality.ConverterJR2Heds;
import de.varylab.discreteconformal.functional.EuclideanNewFunctional;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CAlpha;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CInitialEnergy;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CLambda;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CTheta;
import de.varylab.discreteconformal.unwrapper.numerics.Adapters.CVariable;
import de.varylab.discreteconformal.unwrapper.numerics.MTJDomain;
import de.varylab.discreteconformal.util.UnwrapUtility;

public class EuclideanLayoutTest {

	private CoHDS 	
		hds01 = new CoHDS(),
		hdsCat = new CoHDS();
	private ConverterJR2Heds 
		converter = new ConverterJR2Heds();
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
	public EuclideanNewFunctional<CoVertex, CoEdge, CoFace>
		fun = null;
	private AdapterSet 
		a = AdapterSet.createGenericAdapters();
	
	
	@Before
	public void setUpBeforeClass() throws Exception {
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter());
		
		ReaderOBJ reader01 = new ReaderOBJ();
		Input in01 = new Input("Obj File 01", getClass().getResourceAsStream("tetraflat.obj"));
		SceneGraphComponent c01 = reader01.read(in01);
		IndexedFaceSet ifs01 = (IndexedFaceSet)SceneGraphUtility.getFirstGeometry(c01);
		converter.ifs2heds(ifs01, hds01, a);
		hds01.normalizeCoordinates();
		hdsCat = (CoHDS)HalfedgeIO.readHDS(getClass().getResourceAsStream("cathead.heml"));

		fun = new EuclideanNewFunctional<CoVertex, CoEdge, CoFace>(variable, theta, lambda, alpha, energy);
	}

	@Test
	public void testDoLayout() {
		int n = UnwrapUtility.prepareInvariantDataEuclidean(fun, hds01, a);
		Vector u = new SparseVector(n);
		EuclideanLayout.doLayout(hds01, fun, u);
		for (CoEdge e : hds01.getEdges()) {
			CoVertex s = e.getStartVertex();
			CoVertex t = e.getTargetVertex();
			double l1 = Pn.distanceBetween(s.P, t.P, Pn.EUCLIDEAN);
			double l2 = Pn.distanceBetween(s.T, t.T, Pn.EUCLIDEAN);
			Assert.assertEquals(l1, l2, 1E-11);
		}
	}
	
	
	@Test
	public void testLayout02() throws Exception {
		double[] uDataCatHead = {
			-1.6403228323865509,
			-1.6578647154479722,
			-1.7027277499579634,
			-1.601460767665994,
			-1.6631392800325215,
			-1.6816079958856112,
			-1.667217167609929,
			-1.7117872535525718,
			-1.4902308039019305,
			-1.9239880129680509,
			-2.3393421021580445,
			-2.3677358143759717,
			-2.764693866633633,
			-2.752277788324052,
			-1.568446769965921,
			-2.2714800799942076,
			-2.8942636161289492,
			-2.0593018679664308,
			-2.5393203489584386,
			-2.4908721635477624,
			-2.0190115811624314,
			-1.4206996533916234,
			-1.2470152200202258,
			-2.0189641822231756,
			-1.4035000384813527,
			-1.859887741527121,
			-1.2004511884208808,
			-1.3821594382817353,
			-1.0212662452039216,
			-1.5075394331965881,
			-1.2862854896013192,
			-1.5844814325951575,
			-1.5966593857113534,
			-1.6446045705760637,
			-1.6279061850875252,
			-1.6307280010977314,
			-1.6139788620170148,
			-1.9753629397751886,
			-1.5874246894238542,
			-1.4965326649329456,
			-1.9442050136209825,
			-2.2830039956765473,
			-2.6469279171409195,
			-2.955930761977041,
			-2.7318161277790503,
			-2.3224681485081846,
			-1.9759800892561343,
			-2.328059258861495,
			-2.3120591249606095,
			-1.7160374429888485,
			-2.3717650024161983,
			-1.6903310724188316,
			-1.4653180690613177,
			-1.6854189901696575,
			-1.450109585442535,
			-1.5023579780798617,
			-1.3609286106215883,
			-1.374876635074467,
			-1.4496097245354054,
			-1.5609904790780615,
			-1.5537908593123413,
			-1.4704101403787222,
			-1.2359150948530568,
			-1.318926360273366,
			-1.5651974676583469,
			-1.362212619013763,
			-1.2733205576934323,
			-1.2151983293605852,
			-1.4121378425467326,
			-1.4986019432300326,
			-1.4560338323000843,
			-1.465211486122467,
			-1.0313411542945607,
			-1.2644445584355843,
			-1.2182376839732332,
			-1.4114014853248884,
			-1.6259737748275764,
			-1.4942264779173364,
			-1.4502978903995234,
			-1.5089539832406387,
			-0.997333536262749,
			-1.4360245145328672,
			-1.1863231325591164,
			-1.369712894533045,
			-0.5957729762375399,
			-0.17821262709185479,
			-0.6025863838867218,
			-1.4355885579590946,
			-1.3241330890246312,
			-1.2767799697265272,
			-1.4183135064415322,
			-1.2376111192318393,
			-0.7400916905120952,
			-0.5218975419240646,
			-0.35277381520453077,
			-0.6984654297287117,
			-0.957146262473816,
			-1.0632921582150456,
			-0.8264849780687482,
			-1.2696659621425173,
			-1.2088682813555807,
			-0.7214674525996669,
			-0.5635468912054047,
			-0.8181993194471182,
			-0.29763630140791447,
			-0.21364214633807765,
			-0.12827710497051048,
			-0.3867161891236889,
			-0.32408244679115866,
			-0.9215897885176321,
			-0.5377861837015669,
			-0.23565438869379665,
			-0.9754250713219599,
			-0.6709087875129973,
			-0.37943367850761933,
			-0.20189283637597114,
			-1.205082224804125,
			-0.15131347273232232,
			-1.43756184133558
		};
		Vector uCat = new DenseVector(uDataCatHead);
		MTJDomain dCat = new MTJDomain(uCat);
		// check cat head
		EuclideanLayout.doLayout(hdsCat, fun, uCat);
		for (CoVertex v : hdsCat.getVertices()) {
			if (HalfEdgeUtils.isBoundaryVertex(v)) {
				continue;
			}
			double as = EuclideanLayout.calculateAngleSum(v);
			Assert.assertEquals(2*Math.PI, as, 1E-6);
		}
		for (CoEdge e : hdsCat.getPositiveEdges()) {
			CoVertex s = e.getStartVertex();
			CoVertex t = e.getTargetVertex();
			double tLength = Pn.distanceBetween(s.T, t.T, Pn.EUCLIDEAN);
			double uLength = fun.getNewLength(e, dCat);
			Assert.assertEquals(uLength, tLength, 1E-6);
		}
	}
	
	
}
