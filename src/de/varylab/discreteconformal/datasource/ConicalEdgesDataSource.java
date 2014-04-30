package de.varylab.discreteconformal.datasource;

import de.jreality.geometry.IndexedLineSetUtility;
import de.jreality.scene.Appearance;
import de.jreality.scene.SceneGraphComponent;
import de.jreality.scene.SceneGraphNode;
import de.jreality.shader.CommonAttributes;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.data.DataSourceProvider;
import de.jtem.jrworkspace.plugin.Plugin;
import de.varylab.discreteconformal.adapter.ConeMapAdapter;

public class ConicalEdgesDataSource extends Plugin implements DataSourceProvider {

	private Appearance arcAppearance;


	public ConicalEdgesDataSource() {
		arcAppearance = new Appearance();
		arcAppearance.setAttribute(CommonAttributes.VERTEX_DRAW, false);
		arcAppearance.setAttribute(CommonAttributes.EDGE_DRAW, true);
		arcAppearance.setAttribute(CommonAttributes.FACE_DRAW, false);
	}
	
	
	private class ConicalArcsAdapter extends AbstractAdapter<SceneGraphNode> {

		public ConicalArcsAdapter() {
			super(SceneGraphNode.class, true, false);
		}
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nc) {
			return Edge.class.isAssignableFrom(nc);
		}
		@Override
		public String toString() {
			return "Conical Arcs";
		}
		@Override
		public <
			V extends Vertex<V, E, F>, 
			E extends Edge<V, E, F>, 
			F extends Face<V, E, F>
		> SceneGraphNode getE(E e, AdapterSet a) {
			if (!e.isPositive()) {
				return null;
			}
			SceneGraphComponent root = createArcAtEdge(e, a);
			return root;
		}

		private <
			F extends Face<V, E, F>, 
			E extends Edge<V, E, F>, 
			V extends Vertex<V, E, F>
		> SceneGraphComponent createArcAtEdge(E e, AdapterSet a) {
			SceneGraphComponent c = new SceneGraphComponent(e.toString());
			ConeMapAdapter cma = a.query(ConeMapAdapter.class);
			if(cma == null) {
				return null;
			} else {
				double[][] points =  cma.getArc(e,a,50);
				c.setGeometry(IndexedLineSetUtility.createCurveFromPoints(points, false));
			}
			c.setAppearance(arcAppearance);
			return c;
		}
	}
	
	@Override
	public AdapterSet getDataSources() {
		return new AdapterSet(new ConicalArcsAdapter());
	}

}
