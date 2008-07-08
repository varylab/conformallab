package de.varylab.discreteconformal.frontend.io;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * OFF Reader
 * Written after specification by:
 * http://www.geomview.org/docs/html/geomview_41.html
 */
public final class OFFReader implements IFS {
    private double [][] vertices;
    private double [][] normals;
    private int [][] faces;
    private BufferedReader bufferedReader;
    private String fileType;
    private int vertexCount;
    private int faceCount;
    private boolean hasNormals;
    private int vertexDim;
    private double globalScale;

    /**
     * Constructor takes the filename of an OFF file and instantiously reads it.
     *
     * @param filename
     * @throws Exception
     */
    public OFFReader(String filename) throws Exception {
        this(filename,1);
    }

    /**
     * Constructor takes a stream of an OFF file and instantiously reads it.
     *
     * @param stream
     * @throws Exception
     */
    public OFFReader(InputStream stream) throws Exception {
        this(stream,1);
    }

    /**
     * Constructor takes the filename of an OFF file and instantiously reads it.
     * globalScale scales every vertex by given factor
     *
     * @param filename
     * @param globalScale
     * @throws Exception
     */
    public OFFReader(String filename, double globalScale) throws Exception {
        this.globalScale = globalScale;
        FileInputStream in = new FileInputStream(filename);
        read(in);
    }

    /**
     * Constructor takes a stream of an OFF file and instantiously reads it.
     * globalScale scales every vertex by given factor
     *
     * @param stream
     * @param globalScale
     * @throws Exception
     */
    public OFFReader(InputStream stream, double globalScale) throws Exception {
        this.globalScale = globalScale;
        read(stream);
    }


    public final double[][] getVertices() {
        return vertices;
    }

    public final double[][] getVertexNormals() {
        return normals;
    }

    public final int[][] getFaces() {
        return faces;
    }

    public final int getVertexCount() {
        return vertexCount;
    }

    public final int getFaceCount() {
        return faceCount;
    }

    public final boolean hasVertexNormals() {
        return hasNormals;
    }

    public final int getVertexDim() {
        return vertexDim;
    }

    /**
     * Read the filetype from buffered reader from the header keyword:
     * SPEC:
     * [ST][C][N][4][n]OFF  # Header keyword.
     *
     * @throws Exception
     */
    private void readFileType() throws Exception {
        String readStr = readLine();
        Pattern pattern = Pattern.compile("(.*)OFF");
        Matcher matcher = pattern.matcher(readStr);
        if (!matcher.matches()) {
            bufferedReader.close();
            throw new Exception("File is not in OFF Format");
        }
        fileType = matcher.group(1);
        hasNormals = (fileType.equals("N"));
    }

    /**
     * Read the vertex dimension from file if nOFF format,
     * set to 4 if 4OFF and set to 3 otherwise
     * SPEC:
     * [Ndim]		# Space dimension of vertices, present only if nOFF.
     *
     * @throws Exception
     */
    public final void readVertexDim() throws Exception {
        if (fileType.equals("4"))
            vertexDim = 4;
        else if (fileType.equals("n")) {
            String readStr = readLine();
            vertexDim = Integer.parseInt(readStr);
        } else
            vertexDim = 3;
    }


    /**
     * Read number of vertices, faces and edges
     * SPEC:
     * NVertices  NFaces  NEdges  # NEdges not used or checked.
     *
     * @throws Exception
     */
    private void readIndexCount() throws Exception {
        String readStr = readLine();
        StringTokenizer strToken = new StringTokenizer(readStr);
        vertexCount = Integer.parseInt(strToken.nextToken());
        faceCount = Integer.parseInt(strToken.nextToken());
    }

    /**
     * Read a line from the buffered reader skipping empty lines and truncating comments.
     *
     * @return String of the line read
     * @throws Exception
     */
    private String readLine() throws Exception {
        String readStr;
        readStr = bufferedReader.readLine();
        if (readStr == null) return null;
        readStr = readStr.replaceAll("#.*","").trim();
        if (readStr.equals("")) return readLine();
        return readStr;
    }


    /**
     * Read the vertices (and if available, normals) from the file.
     *
     * @throws Exception
     */
    private void readVertices() throws Exception {
        if (vertexCount == 0) throw new Exception("No vertices in file");
        vertices = new double[vertexCount][3];
        if (hasNormals)
            normals = new double[vertexCount][3];

        StringTokenizer strToken;
        String readStr;
        int i, j;
        for (i = 0; i < vertexCount; i++) {
            readStr = readLine();
            if (readStr == null) throw new EOFException("Expected more vertices (" + vertexCount + ")");
            strToken = new StringTokenizer(readStr);
            if (hasNormals && strToken.countTokens() != vertexDim + 3)
                throw new Exception("Vertex with normals invalid ("+( vertexDim + 3 )+"): " + readStr);
            if (! hasNormals && strToken.countTokens() != vertexDim) throw new Exception("Vertex invalid: " + readStr);
            double[] vertex = new double[vertexDim];
            for (j = 0; j < vertexDim; j++)
                vertex[j] = Float.parseFloat(strToken.nextToken())*globalScale;
            vertices[i] = vertex;
            if (hasNormals) {
                double[] normal = new double[3];
                for (j = 0; j < 3; j++)
                    normal[j] = Float.parseFloat(strToken.nextToken());
                normals[i] = normal;
            }
        }
    }

    /**
     * Read the faces from the file.
     *
     * SPEC:
     * facesN  Vert1 Vert2 ... VertN  [color]
     * (color is currently not supported).
     *
     * @throws Exception
     */
    private void readFaces() throws Exception {
        if (faceCount == 0) return;
        faces = new int[faceCount][];
        StringTokenizer strToken;
        String readStr;
        int faceDim;
        int[] face;
        int i, j;
        for (i = 0; i < faceCount; i++) {
            readStr = readLine();
            if (readStr == null) throw new EOFException("Expected more faces (" + faceCount + ")");
            strToken = new StringTokenizer(readStr);
            faceDim = Integer.parseInt(strToken.nextToken());
            face = new int[faceDim];
            for (j = 0; j < faceDim; j++) {
                face[j] = Integer.parseInt(strToken.nextToken());
            }
            faces[i] = face;
        }
    }

    /**
     * The actual reader takes the stream and parses the data.
     *
     * @param stream
     * @throws Exception
     */
    private void read(InputStream stream) throws Exception {
        faces = null;
        vertices = null;
        normals = null;
        try {
            bufferedReader = new BufferedReader(new InputStreamReader(stream));
            readFileType();
            readVertexDim();
            readIndexCount();
            readVertices();
            readFaces();
            bufferedReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
