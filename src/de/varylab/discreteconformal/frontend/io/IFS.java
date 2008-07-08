package de.varylab.discreteconformal.frontend.io;

public interface IFS {
     public double[][] getVertices();
     public double[][] getVertexNormals();
     public int[][] getFaces();
     public int getVertexCount();
     public int getFaceCount();
     public boolean hasVertexNormals();
     public int getVertexDim();
}
