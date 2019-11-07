<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output method="text"/>

<xsl:template match="/">
  <xsl:apply-templates select="job"/>
</xsl:template>
<xsl:template match="job">
  <xsl:apply-templates select="vars"/>
</xsl:template>
<xsl:template match="vars">
  <xsl:apply-templates select="var"/>
</xsl:template>
<xsl:template match="var[@name='EXPID']">
  <xsl:value-of select="."/>
</xsl:template>

</xsl:stylesheet>
