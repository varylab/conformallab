package de.varylab.discreteconformal.unwrapper;

import java.util.LinkedList;
import java.util.List;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.scene.IndexedFaceSet;
import de.jreality.scene.data.Attribute;
import de.jreality.scene.data.DataList;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.jreality.ConverterHeds2JR;
import de.varylab.discreteconformal.ConformalAdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomEdgeInfo;

public class CircleDomainUnwrapperTest {

	private CoHDS
		hds = null;
	private CoVertex
		v0 = null,
		v1 = null,
		v2 = null,
		v3 = null,
		v4 = null,
		v5 = null,
		v6 = null,
		v7 = null;
	private CoFace
		innerFace01 = null,
		innerFace02 = null;
	private AdapterSet 
		a = new ConformalAdapterSet();
	
	
	@Before
	public void createDataStructure() {
		hds = new CoHDS();
		v0 = hds.addNewVertex();
		v1 = hds.addNewVertex();
		v2 = hds.addNewVertex();
		v3 = hds.addNewVertex();
		v4 = hds.addNewVertex();
		v5 = hds.addNewVertex();
		v6 = hds.addNewVertex();
		v7 = hds.addNewVertex();
		HalfEdgeUtils.constructFaceByVertices(hds, v0, v5, v4);
		HalfEdgeUtils.constructFaceByVertices(hds, v0, v1, v5);
		HalfEdgeUtils.constructFaceByVertices(hds, v1, v6, v5);
		HalfEdgeUtils.constructFaceByVertices(hds, v1, v2, v6);
		HalfEdgeUtils.constructFaceByVertices(hds, v2, v7, v6);
		HalfEdgeUtils.constructFaceByVertices(hds, v2, v3, v7);
		
		innerFace01 = HalfEdgeUtils.constructFaceByVertices(hds, v4, v5, v6).getLeftFace();
		innerFace02 = HalfEdgeUtils.constructFaceByVertices(hds, v6, v7, v4).getLeftFace();
		
		HalfEdgeUtils.constructFaceByVertices(hds, v3, v4, v7);
		HalfEdgeUtils.constructFaceByVertices(hds, v0, v4, v3);
		
		a.set(Position.class, v0, new double[]{0,0});
		 // make it a rectangle
		a.set(Position.class, v1, new double[]{3,0});
		a.set(Position.class, v2, new double[]{5,3});
		a.set(Position.class, v3, new double[]{0,3});
		a.set(Position.class, v4, new double[]{1,1});
		a.set(Position.class, v5, new double[]{2,1});
		a.set(Position.class, v6, new double[]{2,2});
		a.set(Position.class, v7, new double[]{1,2});
	}
	
	
	@Test
	public void testUnwrap() throws Exception {
		CircleDomainUnwrapper unwrapper = new CircleDomainUnwrapper();
		unwrapper.setUnwrapRoot(hds.getVertex(0));
		unwrapper.setGradientTolerance(1E-10);
		unwrapper.setMaxIterations(50);
		unwrapper.unwrap(hds, 0, a);
		
//		TestUtility.display(hds, true);
//		synchronized (this) {
//			wait();	
//		}
//		
		
		double[] vec1 = Rn.subtract(null, v1.T, v0.T);
		double[] vec2 = Rn.subtract(null, v3.T, v0.T);
		double angle1 = Rn.euclideanAngle(vec1, vec2);
		double[] vec3 = Rn.subtract(null, v3.T, v2.T);
		double[] vec4 = Rn.subtract(null, v1.T, v2.T);
		double angle2 = Rn.euclideanAngle(vec3, vec4);
		Assert.assertEquals(Math.PI, angle1 + angle2, 1E-8);
	}
	
	@Test
	public void testUnwrapWithCircularEdges() throws Exception {
		CoEdge ce = HalfEdgeUtils.findEdgeBetweenFaces(innerFace01, innerFace02);
		CustomEdgeInfo info = new CustomEdgeInfo();
		info.circularHoleEdge = true;
		ce.info = info;
		ce.getOppositeEdge().info = info;
		
		CircleDomainUnwrapper unwrapper = new CircleDomainUnwrapper();
		unwrapper.setUnwrapRoot(hds.getVertex(1));
		unwrapper.unwrap(hds, 0, a);
		
		// outer circle
		double[] vec1 = Rn.subtract(null, v1.T, v0.T);
		double[] vec2 = Rn.subtract(null, v3.T, v0.T);
		double angle1 = Rn.euclideanAngle(vec1, vec2);
		double[] vec3 = Rn.subtract(null, v3.T, v2.T);
		double[] vec4 = Rn.subtract(null, v1.T, v2.T);
		double angle2 = Rn.euclideanAngle(vec3, vec4);
		Assert.assertEquals(Math.PI, angle1 + angle2, 1E-8);		
	
		// inner circle
		Rn.subtract(vec1, v5.T, v4.T);
		Rn.subtract(vec2, v7.T, v4.T);
		angle1 = Rn.euclideanAngle(vec1, vec2);
		Rn.subtract(vec3, v7.T, v6.T);
		Rn.subtract(vec4, v5.T, v6.T);
		angle2 = Rn.euclideanAngle(vec3, vec4);
		Assert.assertEquals(Math.PI, angle1 + angle2, 1E-8);	
	}
	
	
	@Test
	public void testStaticUnwrap() throws Exception {
		CoEdge innerEdge = HalfEdgeUtils.findEdgeBetweenFaces(innerFace01, innerFace02);
		TopologyAlgorithms.removeEdge(innerEdge);
		
		int circleIndex = v4.getIndex();
		
		ConverterHeds2JR converter = new ConverterHeds2JR();
		IndexedFaceSet ifs = converter.heds2ifs(hds, a);
		
		int[] zeroBaryFace = {0, 1, 5};
		double[] zeroBaryWeights = {1.0/3, 1.0/3, 1.0/3};		
		List<Integer> circularVerts = new LinkedList<Integer>();
		circularVerts.add(circleIndex);
		CircleDomainUnwrapper.unwrap(ifs, 0, zeroBaryFace, zeroBaryWeights, false, circularVerts);
		DataList tData = ifs.getVertexAttributes(Attribute.TEXTURE_COORDINATES);
		double[][] tArr = tData.toDoubleArrayArray(null);
			
		// outer circle
		double[] vec1 = Rn.subtract(null, tArr[1], tArr[0]);
		double[] vec2 = Rn.subtract(null, tArr[3], tArr[0]);
		double angle1 = Rn.euclideanAngle(vec1, vec2);
		double[] vec3 = Rn.subtract(null, tArr[3], tArr[2]);
		double[] vec4 = Rn.subtract(null, tArr[1], tArr[2]);
		double angle2 = Rn.euclideanAngle(vec3, vec4);
//		Assert.assertEquals(Math.PI, angle1 + angle2, 1E-8);	
		
		Pn.dehomogenize(tArr, tArr);
		// inner circle
		Rn.subtract(vec1, tArr[5], tArr[4]);
		Rn.subtract(vec2, tArr[7], tArr[4]);
		angle1 = Rn.euclideanAngle(vec1, vec2);
		Rn.subtract(vec3, tArr[7], tArr[6]);
		Rn.subtract(vec4, tArr[5], tArr[6]);
		angle2 = Rn.euclideanAngle(vec3, vec4);
		Assert.assertEquals(Math.PI, angle1 + angle2, 1E-8);	
	}
	
	
	public static void main(String[] args) throws Exception {
		new CircleDomainUnwrapperTest().testUnwrap();
	}
	
}
