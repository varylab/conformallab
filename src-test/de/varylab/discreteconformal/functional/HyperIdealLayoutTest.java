package de.varylab.discreteconformal.functional;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;

import java.io.FileOutputStream;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.type.generic.Position4d;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.jtem.jpetsc.Vec;
import de.jtem.jtao.Tao;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.DiscreteMap;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.uniformization.FundamentalPolygon;
import de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility;
import de.varylab.discreteconformal.unwrapper.HyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperIdealApplication;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class HyperIdealLayoutTest {

	private double[]
		ab = {1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 1.7627471360523435, 1.7627471360523428, 1.7627471360523435, 1.7627471360523428, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531};
		
	@BeforeClass
	public static void initPetsc() {
		NativePathUtility.set("native");
		Tao.Initialize();
		LoggingUtility.initLogging();
	}
	
	@Test
	public void testHyperIdealLayout() throws Exception {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiled();
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		Vec u = new Vec(app.getDomainDimension());
		for (int i = 0; i < app.getDomainDimension(); i++) {
			u.setValue(i, ab[i], INSERT_VALUES);
		}
		// initialize angles
		app.evaluateObjectiveAndGradient(u, null);

		Assert.assertEquals(2, HalfEdgeUtils.getGenus(hds));
		
		Set<CoEdge> cutEdges = new LinkedHashSet<>();
		cutEdges.add(hds.getEdge(0));
		cutEdges.add(hds.getEdge(1));
		cutEdges.add(hds.getEdge(3));
		cutEdges.add(hds.getEdge(5));
		cutEdges.add(hds.getEdge(6));
		cutEdges.add(hds.getEdge(8));
		cutEdges.add(hds.getEdge(11));
		cutEdges.add(hds.getEdge(13));
		cutEdges.add(hds.getEdge(14));
		cutEdges.add(hds.getEdge(16));
		cutEdges.add(hds.getEdge(19));
		cutEdges.add(hds.getEdge(21));
		cutEdges.add(hds.getEdge(22));
		cutEdges.add(hds.getEdge(23));
		
		
		CoVertex root = hds.getVertex(3);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
		CuttingUtility.cutAtEdges(cutInfo, cutEdges);
		cutInfo.cutRoot = root;
		Assert.assertEquals(0, HalfEdgeUtils.getGenus(hds));
		
		Map<CoEdge, Double> lMap = new LinkedHashMap<>();
		for (CoEdge e : hds.getPositiveEdges()) {
			if (e.getLeftFace() == null) {
				e = e.getOppositeEdge();
			}
			double l = app.getEdgeLength(e, u);
			lMap.put(e, l);
			lMap.put(e.getOppositeEdge(), l);
		}
		
		// write beta angles to alpha for the layout
		for (CoEdge e : hds.getEdges()) {
			e.setAlpha(e.getBeta());
		}
		
		HyperbolicLayout.doLayout(hds, root, lMap);
		for (CoEdge e : hds.getPositiveEdges()) {
			double[] s = e.getStartVertex().T;
			double[] t = e.getTargetVertex().T;
			double lExpected = lMap.get(e);
			double l = Pn.distanceBetween(s, t, Pn.HYPERBOLIC);
			Assert.assertEquals(lExpected, l, 1E-8);
		}
		
		FundamentalPolygon P = FundamentalPolygonUtility.constructFundamentalPolygon(cutInfo, Pn.HYPERBOLIC);
		System.out.println(P);
		
		ConformalAdapterSet a = new ConformalAdapterSet();
		DiscreteMap map = DataUtility.toDiscreteMap("Lawson Hyper-Ideal Uniformization", hds, a, TexturePosition4d.class, Position4d.class, cutInfo);
		DataIO.writeConformalData(map, new FileOutputStream("test.xml"));
	}
	
}
