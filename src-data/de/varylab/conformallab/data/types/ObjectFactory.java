//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2013.12.09 um 07:42:44 PM CET 
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
    private final static QName _FuchianData_QNAME = new QName("http://www.varylab.com/conformallab/types", "FuchianData");
    private final static QName _ConformalDataList_QNAME = new QName("http://www.varylab.com/conformallab/types", "ConformalDataList");
    private final static QName _DiscreteEmbedding_QNAME = new QName("http://www.varylab.com/conformallab/types", "DiscreteEmbedding");
    private final static QName _HyperEllipticAlgebraicCurve_QNAME = new QName("http://www.varylab.com/conformallab/types", "HyperEllipticAlgebraicCurve");

    /**
     * Create a new ObjectFactory that can be used to create new instances of schema derived classes for package: de.varylab.conformallab.data.types
     * 
     */
    public ObjectFactory() {
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
     * Create an instance of {@link FuchsianData }
     * 
     */
    public FuchsianData createFuchsianData() {
        return new FuchsianData();
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
     * Create an instance of {@link HyperbolicMotion }
     * 
     */
    public HyperbolicMotion createHyperbolicMotion() {
        return new HyperbolicMotion();
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
     * Create an instance of {@link FuchsianGroup }
     * 
     */
    public FuchsianGroup createFuchsianGroup() {
        return new FuchsianGroup();
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
     * Create an instance of {@link FundamentalPolygon }
     * 
     */
    public FundamentalPolygon createFundamentalPolygon() {
        return new FundamentalPolygon();
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
     * Create an instance of {@link JAXBElement }{@code <}{@link FuchsianData }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "FuchianData")
    public JAXBElement<FuchsianData> createFuchianData(FuchsianData value) {
        return new JAXBElement<FuchsianData>(_FuchianData_QNAME, FuchsianData.class, null, value);
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
     * Create an instance of {@link JAXBElement }{@code <}{@link HyperEllipticAlgebraicCurve }{@code >}}
     * 
     */
    @XmlElementDecl(namespace = "http://www.varylab.com/conformallab/types", name = "HyperEllipticAlgebraicCurve")
    public JAXBElement<HyperEllipticAlgebraicCurve> createHyperEllipticAlgebraicCurve(HyperEllipticAlgebraicCurve value) {
        return new JAXBElement<HyperEllipticAlgebraicCurve>(_HyperEllipticAlgebraicCurve_QNAME, HyperEllipticAlgebraicCurve.class, null, value);
    }

}