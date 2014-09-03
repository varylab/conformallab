//
// Diese Datei wurde mit der JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.6 generiert 
// Siehe <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Änderungen an dieser Datei gehen bei einer Neukompilierung des Quellschemas verloren. 
// Generiert: 2014.09.03 um 12:51:40 PM CEST 
//


package de.varylab.conformallab.data.types;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java-Klasse für EmbeddingSelection complex type.
 * 
 * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
 * 
 * <pre>
 * &lt;complexType name="EmbeddingSelection">
 *   &lt;complexContent>
 *     &lt;extension base="{http://www.varylab.com/conformallab/types}ConformalData">
 *       &lt;sequence>
 *         &lt;element name="VertexSelection">
 *           &lt;complexType>
 *             &lt;complexContent>
 *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                 &lt;sequence>
 *                   &lt;element name="Vertex" maxOccurs="unbounded" minOccurs="0">
 *                     &lt;complexType>
 *                       &lt;complexContent>
 *                         &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                           &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *                           &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
 *                         &lt;/restriction>
 *                       &lt;/complexContent>
 *                     &lt;/complexType>
 *                   &lt;/element>
 *                 &lt;/sequence>
 *               &lt;/restriction>
 *             &lt;/complexContent>
 *           &lt;/complexType>
 *         &lt;/element>
 *         &lt;element name="EdgeSelection">
 *           &lt;complexType>
 *             &lt;complexContent>
 *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                 &lt;sequence>
 *                   &lt;element name="Edge" maxOccurs="unbounded" minOccurs="0">
 *                     &lt;complexType>
 *                       &lt;complexContent>
 *                         &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                           &lt;attribute name="face" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *                           &lt;attribute name="vertex1" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *                           &lt;attribute name="vertex2" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *                           &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
 *                         &lt;/restriction>
 *                       &lt;/complexContent>
 *                     &lt;/complexType>
 *                   &lt;/element>
 *                 &lt;/sequence>
 *               &lt;/restriction>
 *             &lt;/complexContent>
 *           &lt;/complexType>
 *         &lt;/element>
 *         &lt;element name="FaceSelection">
 *           &lt;complexType>
 *             &lt;complexContent>
 *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                 &lt;sequence>
 *                   &lt;element name="Face" maxOccurs="unbounded" minOccurs="0">
 *                     &lt;complexType>
 *                       &lt;complexContent>
 *                         &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *                           &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *                           &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
 *                         &lt;/restriction>
 *                       &lt;/complexContent>
 *                     &lt;/complexType>
 *                   &lt;/element>
 *                 &lt;/sequence>
 *               &lt;/restriction>
 *             &lt;/complexContent>
 *           &lt;/complexType>
 *         &lt;/element>
 *       &lt;/sequence>
 *     &lt;/extension>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "EmbeddingSelection", propOrder = {
    "vertexSelection",
    "edgeSelection",
    "faceSelection"
})
public class EmbeddingSelection
    extends ConformalData
{

    @XmlElement(name = "VertexSelection", required = true)
    protected EmbeddingSelection.VertexSelection vertexSelection;
    @XmlElement(name = "EdgeSelection", required = true)
    protected EmbeddingSelection.EdgeSelection edgeSelection;
    @XmlElement(name = "FaceSelection", required = true)
    protected EmbeddingSelection.FaceSelection faceSelection;

    /**
     * Ruft den Wert der vertexSelection-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link EmbeddingSelection.VertexSelection }
     *     
     */
    public EmbeddingSelection.VertexSelection getVertexSelection() {
        return vertexSelection;
    }

    /**
     * Legt den Wert der vertexSelection-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link EmbeddingSelection.VertexSelection }
     *     
     */
    public void setVertexSelection(EmbeddingSelection.VertexSelection value) {
        this.vertexSelection = value;
    }

    /**
     * Ruft den Wert der edgeSelection-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link EmbeddingSelection.EdgeSelection }
     *     
     */
    public EmbeddingSelection.EdgeSelection getEdgeSelection() {
        return edgeSelection;
    }

    /**
     * Legt den Wert der edgeSelection-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link EmbeddingSelection.EdgeSelection }
     *     
     */
    public void setEdgeSelection(EmbeddingSelection.EdgeSelection value) {
        this.edgeSelection = value;
    }

    /**
     * Ruft den Wert der faceSelection-Eigenschaft ab.
     * 
     * @return
     *     possible object is
     *     {@link EmbeddingSelection.FaceSelection }
     *     
     */
    public EmbeddingSelection.FaceSelection getFaceSelection() {
        return faceSelection;
    }

    /**
     * Legt den Wert der faceSelection-Eigenschaft fest.
     * 
     * @param value
     *     allowed object is
     *     {@link EmbeddingSelection.FaceSelection }
     *     
     */
    public void setFaceSelection(EmbeddingSelection.FaceSelection value) {
        this.faceSelection = value;
    }


    /**
     * <p>Java-Klasse für anonymous complex type.
     * 
     * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
     * 
     * <pre>
     * &lt;complexType>
     *   &lt;complexContent>
     *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *       &lt;sequence>
     *         &lt;element name="Edge" maxOccurs="unbounded" minOccurs="0">
     *           &lt;complexType>
     *             &lt;complexContent>
     *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *                 &lt;attribute name="face" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
     *                 &lt;attribute name="vertex1" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
     *                 &lt;attribute name="vertex2" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
     *                 &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
     *               &lt;/restriction>
     *             &lt;/complexContent>
     *           &lt;/complexType>
     *         &lt;/element>
     *       &lt;/sequence>
     *     &lt;/restriction>
     *   &lt;/complexContent>
     * &lt;/complexType>
     * </pre>
     * 
     * 
     */
    @XmlAccessorType(XmlAccessType.FIELD)
    @XmlType(name = "", propOrder = {
        "edges"
    })
    public static class EdgeSelection {

        @XmlElement(name = "Edge")
        protected List<EmbeddingSelection.EdgeSelection.Edge> edges;

        /**
         * Gets the value of the edges property.
         * 
         * <p>
         * This accessor method returns a reference to the live list,
         * not a snapshot. Therefore any modification you make to the
         * returned list will be present inside the JAXB object.
         * This is why there is not a <CODE>set</CODE> method for the edges property.
         * 
         * <p>
         * For example, to add a new item, do as follows:
         * <pre>
         *    getEdges().add(newItem);
         * </pre>
         * 
         * 
         * <p>
         * Objects of the following type(s) are allowed in the list
         * {@link EmbeddingSelection.EdgeSelection.Edge }
         * 
         * 
         */
        public List<EmbeddingSelection.EdgeSelection.Edge> getEdges() {
            if (edges == null) {
                edges = new ArrayList<EmbeddingSelection.EdgeSelection.Edge>();
            }
            return this.edges;
        }


        /**
         * <p>Java-Klasse für anonymous complex type.
         * 
         * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
         * 
         * <pre>
         * &lt;complexType>
         *   &lt;complexContent>
         *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
         *       &lt;attribute name="face" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
         *       &lt;attribute name="vertex1" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
         *       &lt;attribute name="vertex2" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
         *       &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
         *     &lt;/restriction>
         *   &lt;/complexContent>
         * &lt;/complexType>
         * </pre>
         * 
         * 
         */
        @XmlAccessorType(XmlAccessType.FIELD)
        @XmlType(name = "")
        public static class Edge {

            @XmlAttribute(name = "face", required = true)
            protected int face;
            @XmlAttribute(name = "vertex1", required = true)
            protected int vertex1;
            @XmlAttribute(name = "vertex2", required = true)
            protected int vertex2;
            @XmlAttribute(name = "channel")
            protected Integer channel;

            /**
             * Ruft den Wert der face-Eigenschaft ab.
             * 
             */
            public int getFace() {
                return face;
            }

            /**
             * Legt den Wert der face-Eigenschaft fest.
             * 
             */
            public void setFace(int value) {
                this.face = value;
            }

            /**
             * Ruft den Wert der vertex1-Eigenschaft ab.
             * 
             */
            public int getVertex1() {
                return vertex1;
            }

            /**
             * Legt den Wert der vertex1-Eigenschaft fest.
             * 
             */
            public void setVertex1(int value) {
                this.vertex1 = value;
            }

            /**
             * Ruft den Wert der vertex2-Eigenschaft ab.
             * 
             */
            public int getVertex2() {
                return vertex2;
            }

            /**
             * Legt den Wert der vertex2-Eigenschaft fest.
             * 
             */
            public void setVertex2(int value) {
                this.vertex2 = value;
            }

            /**
             * Ruft den Wert der channel-Eigenschaft ab.
             * 
             * @return
             *     possible object is
             *     {@link Integer }
             *     
             */
            public int getChannel() {
                if (channel == null) {
                    return  0;
                } else {
                    return channel;
                }
            }

            /**
             * Legt den Wert der channel-Eigenschaft fest.
             * 
             * @param value
             *     allowed object is
             *     {@link Integer }
             *     
             */
            public void setChannel(Integer value) {
                this.channel = value;
            }

        }

    }


    /**
     * <p>Java-Klasse für anonymous complex type.
     * 
     * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
     * 
     * <pre>
     * &lt;complexType>
     *   &lt;complexContent>
     *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *       &lt;sequence>
     *         &lt;element name="Face" maxOccurs="unbounded" minOccurs="0">
     *           &lt;complexType>
     *             &lt;complexContent>
     *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *                 &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
     *                 &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
     *               &lt;/restriction>
     *             &lt;/complexContent>
     *           &lt;/complexType>
     *         &lt;/element>
     *       &lt;/sequence>
     *     &lt;/restriction>
     *   &lt;/complexContent>
     * &lt;/complexType>
     * </pre>
     * 
     * 
     */
    @XmlAccessorType(XmlAccessType.FIELD)
    @XmlType(name = "", propOrder = {
        "faces"
    })
    public static class FaceSelection {

        @XmlElement(name = "Face")
        protected List<EmbeddingSelection.FaceSelection.Face> faces;

        /**
         * Gets the value of the faces property.
         * 
         * <p>
         * This accessor method returns a reference to the live list,
         * not a snapshot. Therefore any modification you make to the
         * returned list will be present inside the JAXB object.
         * This is why there is not a <CODE>set</CODE> method for the faces property.
         * 
         * <p>
         * For example, to add a new item, do as follows:
         * <pre>
         *    getFaces().add(newItem);
         * </pre>
         * 
         * 
         * <p>
         * Objects of the following type(s) are allowed in the list
         * {@link EmbeddingSelection.FaceSelection.Face }
         * 
         * 
         */
        public List<EmbeddingSelection.FaceSelection.Face> getFaces() {
            if (faces == null) {
                faces = new ArrayList<EmbeddingSelection.FaceSelection.Face>();
            }
            return this.faces;
        }


        /**
         * <p>Java-Klasse für anonymous complex type.
         * 
         * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
         * 
         * <pre>
         * &lt;complexType>
         *   &lt;complexContent>
         *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
         *       &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
         *       &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
         *     &lt;/restriction>
         *   &lt;/complexContent>
         * &lt;/complexType>
         * </pre>
         * 
         * 
         */
        @XmlAccessorType(XmlAccessType.FIELD)
        @XmlType(name = "")
        public static class Face {

            @XmlAttribute(name = "index", required = true)
            protected int index;
            @XmlAttribute(name = "channel")
            protected Integer channel;

            /**
             * Ruft den Wert der index-Eigenschaft ab.
             * 
             */
            public int getIndex() {
                return index;
            }

            /**
             * Legt den Wert der index-Eigenschaft fest.
             * 
             */
            public void setIndex(int value) {
                this.index = value;
            }

            /**
             * Ruft den Wert der channel-Eigenschaft ab.
             * 
             * @return
             *     possible object is
             *     {@link Integer }
             *     
             */
            public int getChannel() {
                if (channel == null) {
                    return  0;
                } else {
                    return channel;
                }
            }

            /**
             * Legt den Wert der channel-Eigenschaft fest.
             * 
             * @param value
             *     allowed object is
             *     {@link Integer }
             *     
             */
            public void setChannel(Integer value) {
                this.channel = value;
            }

        }

    }


    /**
     * <p>Java-Klasse für anonymous complex type.
     * 
     * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
     * 
     * <pre>
     * &lt;complexType>
     *   &lt;complexContent>
     *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *       &lt;sequence>
     *         &lt;element name="Vertex" maxOccurs="unbounded" minOccurs="0">
     *           &lt;complexType>
     *             &lt;complexContent>
     *               &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
     *                 &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
     *                 &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
     *               &lt;/restriction>
     *             &lt;/complexContent>
     *           &lt;/complexType>
     *         &lt;/element>
     *       &lt;/sequence>
     *     &lt;/restriction>
     *   &lt;/complexContent>
     * &lt;/complexType>
     * </pre>
     * 
     * 
     */
    @XmlAccessorType(XmlAccessType.FIELD)
    @XmlType(name = "", propOrder = {
        "vertices"
    })
    public static class VertexSelection {

        @XmlElement(name = "Vertex")
        protected List<EmbeddingSelection.VertexSelection.Vertex> vertices;

        /**
         * Gets the value of the vertices property.
         * 
         * <p>
         * This accessor method returns a reference to the live list,
         * not a snapshot. Therefore any modification you make to the
         * returned list will be present inside the JAXB object.
         * This is why there is not a <CODE>set</CODE> method for the vertices property.
         * 
         * <p>
         * For example, to add a new item, do as follows:
         * <pre>
         *    getVertices().add(newItem);
         * </pre>
         * 
         * 
         * <p>
         * Objects of the following type(s) are allowed in the list
         * {@link EmbeddingSelection.VertexSelection.Vertex }
         * 
         * 
         */
        public List<EmbeddingSelection.VertexSelection.Vertex> getVertices() {
            if (vertices == null) {
                vertices = new ArrayList<EmbeddingSelection.VertexSelection.Vertex>();
            }
            return this.vertices;
        }


        /**
         * <p>Java-Klasse für anonymous complex type.
         * 
         * <p>Das folgende Schemafragment gibt den erwarteten Content an, der in dieser Klasse enthalten ist.
         * 
         * <pre>
         * &lt;complexType>
         *   &lt;complexContent>
         *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
         *       &lt;attribute name="index" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
         *       &lt;attribute name="channel" type="{http://www.w3.org/2001/XMLSchema}int" default="0" />
         *     &lt;/restriction>
         *   &lt;/complexContent>
         * &lt;/complexType>
         * </pre>
         * 
         * 
         */
        @XmlAccessorType(XmlAccessType.FIELD)
        @XmlType(name = "")
        public static class Vertex {

            @XmlAttribute(name = "index", required = true)
            protected int index;
            @XmlAttribute(name = "channel")
            protected Integer channel;

            /**
             * Ruft den Wert der index-Eigenschaft ab.
             * 
             */
            public int getIndex() {
                return index;
            }

            /**
             * Legt den Wert der index-Eigenschaft fest.
             * 
             */
            public void setIndex(int value) {
                this.index = value;
            }

            /**
             * Ruft den Wert der channel-Eigenschaft ab.
             * 
             * @return
             *     possible object is
             *     {@link Integer }
             *     
             */
            public int getChannel() {
                if (channel == null) {
                    return  0;
                } else {
                    return channel;
                }
            }

            /**
             * Legt den Wert der channel-Eigenschaft fest.
             * 
             * @param value
             *     allowed object is
             *     {@link Integer }
             *     
             */
            public void setChannel(Integer value) {
                this.channel = value;
            }

        }

    }

}
