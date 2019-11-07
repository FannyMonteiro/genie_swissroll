
JC=/usr/bin/javac
JAR=/usr/bin/jar
JAVA=/usr/bin/java
CLASSPATH=/usr/share/java/

%.class : %.java
	$(JC) $*.java

# The first entry is the main class
MANIFEST=manifest.mf
MAIN_CLASS=GenieGUI
SOURCEFILES = $(MAIN_CLASS).java \
	AdjstblParam.java \
	Colours.java \
	DefnParam.java \
	GenieFileReader.java \
	GlobePlotter.java \
	GraphPlotter.java \
	HTMLFiles.java \
	JobParam.java \
	KButton.java \
	KSlider.java \
	MathTest.java \
	ParamSetting.java \
	TextDraw.java \
	TextDrawPanel.java \
	TimeSeriesFrame.java \
	Utils.java \
	XMLReader.java \
	XSLT.java \
	util/misc/CausedIOException.java \
	util/xml/XMLException.java \
	util/xml/XML.java \
	util/xml/XSL.java
CLASSES = $(SOURCEFILES:.java=.class)
ALLCLASSES = $(SOURCEFILES:.java=*.class)

JARFILE=$(MAIN_CLASS).jar

all: $(CLASSES) $(MANIFEST)
	$(JAR) cvmf $(MANIFEST) $(JARFILE) $(ALLCLASSES)

.PHONY: clean
clean:
	rm -f $(ALLCLASSES) $(JARFILE)

# dependencies

$(MANIFEST): makefile
	echo "Main-Class: "$(MAIN_CLASS) > $(MANIFEST)
	echo "Class-Path: "$(CLASSPATH)"/serializer.jar "$(CLASSPATH)"/xml-apis.jar "$(CLASSPATH)"/xercesImpl.jar "$(CLASSPATH)"/xalan.jar" >> $(MANIFEST)

test: all
	echo
	echo "JC="$(JC)
	echo "JAVA="$(JAVA)
	echo "JC & JAVA vals may need changing see makefile"
	echo
	cp $(JARFILE) ../../test
	cd ../../test; $(JAVA) -jar $(JARFILE)

ENGAGEtest: all
	echo
	echo "JC="$(JC)
	echo "JAVA="$(JAVA)
	echo "JC & JAVA vals may need changing see makefile"
	echo
	./ENGAGEtest.create_runData.xml.sh
	cp $(JARFILE) ../../test
	cd ../../test; $(JAVA) -jar $(JARFILE)

