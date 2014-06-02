//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.05.30 um 12:40:20 PM CEST 
//


package de.varylab.conformallab.data.types;

import javax.xml.bind.JAXBElement;
import javax.xml.bind.annotation.XmlElementDecl;
import javax.xml.bind.annotation.XmlRegistry;
import javax.xml.namespace.QName;


/**
 * This object contains factory methods for each 
 * Java content interface and Java element interface 
 * generated in the de.varylab.conformallab.data.types package. 
 * <p>An ObjectFactory allows you to programatically 
 * construct new instances of the Java representation 
 * for XML content. The Java representation of XML 
 * content can consist of schema derived interfaces 
 * and classes representing the binding of schema 
 * type definitions, element declarations and model 
 * groups.  Factory methods for each of these are 
 * provided in this class.
 * 
 */
@XmlRegistry
public class ObjectFactory {

    private final static QName _DiscreteMetric_QNAME = new QName("http://www.varylab.com/conformallab/types", "DiscreteMetric");
    private final static QName _SchottkyData_QNAME = new QName("http://www.varylab.com/conformallab/types", "SchottkyData");
    private final static QName _DiscreteMap_QNAME = new QName("http://www.varylab.com/conformallab/types", "DiscreteMap");
    private final static QName _UniformizationData_QNAME = new QName("http://www.varylab.com/conformallab/types", "UniformizationData");
    private final static QName _ConformalDataList_QNAME = new QName("http://www.varylab.com/conformallab/types", "ConformalDataList");
    private final static QName _DiscreteEmbedding_QNAME = new QName("http://www.varylab.com/conformallab/types", "DiscreteEmbedding");
    private final static QName _EmbeddingSelection_QNAME = new QName("http://www.varylab.com/conformallab/types", "EmbeddingSelection");
    private final static QName _HyperEllipticAlgebraicCurve_QNAME = new QName("http://www.varylab.com/conformallab/types", "HyperEllipticAlgebraicCurve");

    /**
     * Create a new ObjectFactory that can be used to create new instances of schema derived classes for package: de.varylab.conformallab.data.types
     * 
     */
    public ObjectFactory() {
    }

    /**
     * Create an instance of {@link EmbeddingSelection }
     * 
     */
    public EmbeddingSelection createEmbeddingSelection() {
        return new EmbeddingSelection();
    }

    /**
     * Create an instance of {@link EmbeddingSelection.FaceSelection }
     * 
     */
    public EmbeddingSelection.FaceSelection createEmbeddingSelectionFaceSelection() {
        return new EmbeddingSelection.FaceSelection();
    }

    /**
     * Create an instance of {@link EmbeddingSelection.EdgeSelection }
     * 
     */
    public EmbeddingSelection.EdgeSelection createEmbeddingSelectionEdgeSelection() {
        return new EmbeddingSelection.EdgeSelection();
    }

    /**
     * Create an instance of {@link EmbeddingSelection.VertexSelection }
     * 
     */
    public EmbeddingSelection.VertexSelection createEmbeddingSelectionVertexSelection() {
        return new EmbeddingSelection.VertexSelection();
    }

    /**
     * Create an instance of {@link SchottkyData }
     * 
     */
    public SchottkyData createSchottkyData() {
        return new SchottkyData();
    }

    /**
     * Create an instance of {@link DiscreteMetric }
     * 
     */
    public DiscreteMetric createDiscreteMetric() {
        return new DiscreteMetric();
    }

    /**
     * Create an instance of {@link UniformizationData }
     * 
     */
    public UniformizationData createUniformizationData() {
        return new UniformizationData();
    }

    /**
     * Create an instance of {@link DiscreteMap }
     * 
     */
    public DiscreteMap createDiscreteMap() {
        return new DiscreteMap();
    }

    /**
     * Create an instance of {@link DiscreteEmbedding }
     * 
     */
    public DiscreteEmbedding createDiscreteEmbedding() {
        return new DiscreteEmbedding();
    }

    /**
     * Create an instance of {@link ConformalDataList }
     * 
     */
    public ConformalDataList createConformalDataList() {
        return new ConformalDataList();
    }

    /**
     * Create an instance of {@link HyperEllipticAlgebraicCurve }
     * 
     */
    public HyperEllipticAlgebraicCurve createHyperEllipticAlgebraicCurve() {
        return new HyperEllipticAlgebraicCurve();
    }

    /**
     * Create an instance of {@link VertexIdentification }
     * 
     */
    public VertexIdentification createVertexIdentification() {
        return new VertexIdentification();
    }

    /**
     * Create an instance of {@link EmbeddedTriangle }
     * 
     */
    public EmbeddedTriangle createEmbeddedTriangle() {
        return new EmbeddedTriangle();
    }

    /**
     * Create an instance of {@link SchottkyGenerator }
     * 
     */
    public SchottkyGenerator createSchottkyGenerator() {
        return new SchottkyGenerator();
    }

    /**
     * Create an instance of {@link Complex }
     * 
     */
    public Complex createComplex() {
        return new Complex();
    }

    /**
     * Create an instance of {@link IsometryPSL2R }
     * 
     */
    public IsometryPSL2R createIsometryPSL2R() {
        return new IsometryPSL2R();
    }

    /**
     * Create an instance of {@link EmbeddedVertex }
     * 
     */
    public EmbeddedVertex createEmbeddedVertex() {
        return new EmbeddedVertex();
    }

    /**
     * Create an instance of {@link FundamentalEdge }
     * 
     */
    public FundamentalEdge createFundamentalEdge() {
        return new FundamentalEdge();
    }

    /**
     * Create an instance of {@link MetricEdge }
     * 
     */
    public MetricEdge createMetricEdge() {
        return new MetricEdge();
    }

    /**
     * Create an instance of {@link Moebius }
     * 
     */
    public Moebius createMoebius() {
        return new Moebius();
    }

    /**
     * Create an instance of {@link ConformalData }
     * 
     */
    public ConformalData createConformalData() {
        return new ConformalData();
    }

    /**
     * Create an instance of {@link MetricTriangle }
     * 
     */
    public MetricTriangle createMetricTriangle() {
        return new MetricTriangle();
    }

    /**
     * Create an instance of {@link Circle }
     * 
     */
    public Circle createCircle() {
        return new Circle();
    }

    /**
     * Create an instance of {@link FundamentalVertex }
     * 
     */
    public FundamentalVertex createFundamentalVertex() {
        return new FundamentalVertex();
    }

    /**
     * Create an instance of {@link FundamentalPolygon }
     * 
     */
    public FundamentalPolygon createFundamentalPolygon() {
        return new FundamentalPolygon();
    }

    /**
     * Create an instance of {@link UniformizingGroup }
     * 
     */
    public UniformizingGroup createUniformizingGroup() {
        return new UniformizingGroup();
    }

    /**
     * Create an instance of {@link EmbeddingSelection.FaceSelection.Face }
     * 
     */
    public EmbeddingSelection.FaceSelection.Face createEmbeddingSelectionFaceSelectionFace() {
        return new EmbeddingSelection.FaceSelection.Face();
    }

    /**
     * Create an instance of {@link EmbeddingSelection.EdgeSelection.Edge }
     * 
     */
    public EmbeddingSelection.EdgeSelection.Edge createEmbeddingSelectionEdgeSelectionEdge() {
        return new EmbeddingSelection.EdgeSelection.Edge();
    }

    /**
     * Create an instance of {@link EmbeddingSelection.VertexSelection.Vertex }
     * 
     */
    public EmbeddingSelection.VertexSelection.Vertex createEmbeddingSelectionVertexSelectionVertex() {
        return new EmbeddingSelection.VertexSelection.Vertex();
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link DiscreteMetric }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "DiscreteMetric")
    public JAXBElement<DiscreteMetric> createDiscreteMetric(DiscreteMetric value) {
        return new JAXBElement<DiscreteMetric>(_DiscreteMetric_QNAME, DiscreteMetric.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link SchottkyData }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "SchottkyData")
    public JAXBElement<SchottkyData> createSchottkyData(SchottkyData value) {
        return new JAXBElement<SchottkyData>(_SchottkyData_QNAME, SchottkyData.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link DiscreteMap }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "DiscreteMap")
    public JAXBElement<DiscreteMap> createDiscreteMap(DiscreteMap value) {
        return new JAXBElement<DiscreteMap>(_DiscreteMap_QNAME, DiscreteMap.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link UniformizationData }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "UniformizationData")
    public JAXBElement<UniformizationData> createUniformizationData(UniformizationData value) {
        return new JAXBElement<UniformizationData>(_UniformizationData_QNAME, UniformizationData.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link ConformalDataList }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "ConformalDataList")
    public JAXBElement<ConformalDataList> createConformalDataList(ConformalDataList value) {
        return new JAXBElement<ConformalDataList>(_ConformalDataList_QNAME, ConformalDataList.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link DiscreteEmbedding }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "DiscreteEmbedding")
    public JAXBElement<DiscreteEmbedding> createDiscreteEmbedding(DiscreteEmbedding value) {
        return new JAXBElement<DiscreteEmbedding>(_DiscreteEmbedding_QNAME, DiscreteEmbedding.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link EmbeddingSelection }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "EmbeddingSelection")
    public JAXBElement<EmbeddingSelection> createEmbeddingSelection(EmbeddingSelection value) {
        return new JAXBElement<EmbeddingSelection>(_EmbeddingSelection_QNAME, EmbeddingSelection.class, null, value);
    }

    /**
     * Create an instance of {@link JAXBElement }{@code <}{@link HyperEllipticAlgebraicCurve }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "HyperEllipticAlgebraicCurve")
    public JAXBElement<HyperEllipticAlgebraicCurve> createHyperEllipticAlgebraicCurve(HyperEllipticAlgebraicCurve value) {
        return new JAXBElement<HyperEllipticAlgebraicCurve>(_HyperEllipticAlgebraicCurve_QNAME, HyperEllipticAlgebraicCurve.class, null, value);
    }

}
