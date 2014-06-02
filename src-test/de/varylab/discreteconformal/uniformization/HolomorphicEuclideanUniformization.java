package de.varylab.discreteconformal.uniformization;

import static de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility.constructFundamentalPolygon;

import java.io.InputStream;
import java.util.Iterator;
import java.util.logging.Logger;

import org.junit.BeforeClass;
import org.junit.Test;

import cern.colt.Arrays;
import de.jreality.math.FactoredMatrix;
import de.jreality.math.Pn;
import de.jreality.util.NativePathUtility;
import de.jtem.discretegroup.core.DiscreteGroup;
import de.jtem.discretegroup.core.DiscreteGroupElement;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.selection.Selection;
import de.varylab.conformallab.data.DataIO;
import de.varylab.conformallab.data.DataUtility;
import de.varylab.conformallab.data.types.ConformalDataList;
import de.varylab.conformallab.data.types.DiscreteEmbedding;
import de.varylab.conformallab.data.types.DiscreteMap;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomVertexInfo;
import de.varylab.discreteconformal.logging.LoggingUtility;
import de.varylab.discreteconformal.unwrapper.EuclideanUnwrapperPETSc;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class HolomorphicEuclideanUniformization {

	private Logger
		log = Logger.getLogger(HolomorphicEuclideanUniformization.class.getName());
	
	@BeforeClass
	public static void initPetsc() {
		NativePathUtility.set("native");
		LoggingUtility.initLogging();
	}
	
	@Test
	public void testCreateGroup() throws Exception {
		InputStream in = getClass().getResourceAsStream("HolomorphicEuclideanUniformization_blueprint.xml");
		ConformalDataList list = DataIO.readConformalDataList(in);
		DiscreteMap map = (DiscreteMap) list.getData().get(1);
		CuttingInfo<CoVertex, CoEdge, CoFace> cuts = new CuttingInfo<>();
		CoHDS hds = DataUtility.toHDS(map.getDomain(), cuts);
		for (CoVertex v : hds.getVertices()) {
			v.T = v.P.clone();
		}
		FundamentalPolygon P = constructFundamentalPolygon(cuts, Pn.EUCLIDEAN);
		log.info(P.toString());
		DiscreteGroup G = P.getDiscreteGroup();
		for (DiscreteGroupElement g : G.getGenerators()) {
			FactoredMatrix fm = new FactoredMatrix(g.getArray());
			log.info(fm.getRotation().toString());
			double sin = fm.getEntry(0, 1);
			log.info("Angle: " + Math.toDegrees(Math.asin(sin)));
		}
	}

	@Test
	public void testUniformizeEuclidean() throws Exception {
		InputStream in = getClass().getResourceAsStream("HolomorphicEuclideanUniformization_model.xml");
		DiscreteEmbedding de = (DiscreteEmbedding)DataIO.readConformalData(in);
		CuttingInfo<CoVertex, CoEdge, CoFace> cuts = new CuttingInfo<>();
		CoHDS hds = DataUtility.toHDS(de, cuts);
		Selection s = DataUtility.toSelection(de.getSelection(), hds);
		log.info("HDS: " + hds);
		log.info("Selected Vertices: " + s.getVertices());
		log.info("Selected Edges: " + s.getEdges().size());
		Iterator<CoVertex> vIterator = s.getVertices(hds).iterator();
		CoVertex v0 = vIterator.next();
		CoVertex v1 = vIterator.next();
		v0.info = new CustomVertexInfo();
		v1.info = new CustomVertexInfo();
		v0.info.useCustomTheta = true;
		v1.info.useCustomTheta = true;
		v0.info.theta = 4 * Math.PI;
		v1.info.theta = 4 * Math.PI;
		
		// tweak singulatiry position
		log.info("singularity pos " + v0 + ": " + Arrays.toString(v0.P));
		v0.P[1] += 0.033;
		
		// uniformize
		EuclideanUnwrapperPETSc unwrapper = new EuclideanUnwrapperPETSc();
		unwrapper.setCutGraph(s.getEdges(hds));
		AdapterSet aSet = new ConformalAdapterSet();
		unwrapper.unwrap(hds, 2, aSet);
		
		// group
		cuts = unwrapper.getCutInfo();
		FundamentalPolygon P = constructFundamentalPolygon(cuts, Pn.EUCLIDEAN);
		DiscreteGroup G = P.getDiscreteGroup();
		for (DiscreteGroupElement g : G.getGenerators()) {
			FactoredMatrix fm = new FactoredMatrix(g.getArray());
			log.info(fm.getRotation().toString());
			double sin = fm.getEntry(0, 1);
			log.info("Angle: " + Math.toDegrees(Math.asin(sin)));
		}
	}
	
	
}
