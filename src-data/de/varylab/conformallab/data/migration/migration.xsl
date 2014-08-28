<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet 
	version="1.0" 
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns:vd="http://www.varylab.com/conformallab/types"
>
	<xsl:output indent="no" method="xml" omit-xml-declaration="no"/>
	
	<!-- default copy template -->
	<xsl:template match="*">
		<xsl:copy-of select="."/>
	</xsl:template>
	
	<!-- to version 1 -->
	
	<xsl:template match="vd:DiscreteMap[not(@version)]">
		<xsl:copy>
			<xsl:copy-of select="@*"/>
			<xsl:attribute name="version">1</xsl:attribute>
			<xsl:apply-templates select="child::node()"/>
		</xsl:copy>
	</xsl:template>
	
	<xsl:template match="vd:Domain[not(@version)] | vd:Image[not(@version)]">
		<xsl:copy>
			<xsl:copy-of select="@*"/>
			<xsl:attribute name="version">1</xsl:attribute>
			<xsl:apply-templates select="child::node()"/>
		</xsl:copy>
	</xsl:template>
	

	
</xsl:stylesheet>