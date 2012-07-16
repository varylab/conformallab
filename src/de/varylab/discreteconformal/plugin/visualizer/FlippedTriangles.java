package de.varylab.discreteconformal.plugin.visualizer;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureFieldMin;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.halfedgetools.plugin.data.DataSourceProvider;
import de.jtem.jrworkspace.plugin.Plugin;
import de.varylab.discreteconformal.adapter.TriangleOrientationAdapter;
import de.varylab.discreteconformal.unwrapper.isothermic.IsothermicUtility;

public class FlippedTriangles extends Plugin implements DataSourceProvider {

	private TriangleOrientationAdapter adapter = new TriangleOrientationAdapter();	
	
	@Override
	public AdapterSet getDataSources() {
		return new AdapterSet(adapter);
	}

}
