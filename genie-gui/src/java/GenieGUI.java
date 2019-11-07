import java.awt.*;
import java.applet.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
import java.io.*;
import java.util.*;

/* ////////

Component name conventions:

nameB.  KButton (KButton.java)
nameCB. JCheckBox
nameChoice. Choice
nameTF. TextField
nameS. Scrollbar
nameMI. MenuItem
nameTDP. TextDrawPanel (TextDrawPanel.java)
nameMenu. Menu

/////////*/

public class GenieGUI extends Applet implements Runnable{
  boolean isStandalone=false, first=true;

// GENIE parameters/variables

///////////////////////////
	public static int iMax=36, jMax=36, iMax1=iMax-1, jMax1=jMax-1;

///////////////////////
// Arrays which may need adjusting:
//   resize(w,h) will need calling
//   GlobePlotter.java calc function will need modifying

	double timeVal[]=new double[200000] // Increase size for long runs
    , yVals[]=new double[200000]
		, graphVal[]=new double[200000]
   ;
	
	public static double mapVal[][]=new double[iMax][jMax];
	public static double initialMapVal[][]=new double[iMax][jMax];
	
	public static boolean land[][]=new boolean[iMax][jMax];
	public static boolean snow[][]=new boolean[iMax][jMax];
	public static boolean seaBed[][][]=new boolean[iMax][jMax][8];

	Color mapColrs[][]=new Color[iMax][jMax];
// Coastline points may need increasing with map size, 
//   GlobePlotter.java will also need changing
	public int[] clStartX=new int[1000], clStartY=new int[1000], clEndX=new int[1000], clEndY=new int[1000];

//////////////////////////
	
	static String genieOutputDir, genieRunDir, genieInitialDir
    , genieCommand
		, timeSFile, mapFile, initialMapFile;
	public static int colNum, nHeaderLines, nCols;
	
	static double time, compareTime;
	int nReadings=0;
	int nYVals=0;
	int xScale=1;
	int runLength=1000;
	
	public static double mapMin=1e100, mapMax=-1e100;
	double fireStart, fireEnd;
	
/* Variables used in OU gui
	double depths[]={5e3,3.58e3,2.52e3,1.74e3,1.16e3,0.729e3,0.411e3,0.175e3};
	double windU[][]=new double[iMax][jMax];
	double windV[][]=new double[iMax][jMax];
	double viMin[]=new double[27];
	double viMax[]=new double[27];
	
	double velXVal[][][]=new double[iMax][jMax][8];
	double vXMin[]=new double[9], vXMax[]=new double[9];
	double velYVal[][][]=new double[iMax][jMax][8];
	double vYMin[]=new double[9], vYMax[]=new double[9];
	double velScale=15;
//*/
	
	///// Coastline variables
	public int nCoastlines;
	
	///// Globe plotting variables
	double psi=0, psi0=psi, cp=1, sp=0;
	double eps=Math.PI, eps0=eps, ce=1, se=0;
	int Y0=0, X0=0;

	static final double M_PI=Math.PI;

	double grad;
	
	/////// Global variables to control graph dimensions
	
	double xMin, xMax, yMin, yMax;
	
	/////////////////////////////////
	
	int mapDPs=2, graphDPs=2;

	// HTML file variables
	String folderName, HTMLlinks;
	String tableHeader[]={" ", " "};
	int oldAvIdx, oldSeasIdx;
	char qt='"';

///////////////////////////

//   Interface varaiables/components

	static int minWidth=800, minHeight=600;

	static KButton  startB, pauseB, stopB
		, paramsB, paramsOKB, paramsCanB, paramsDefB
		, emsB
		, exptDirB, cOKB, cCanB
		, changeB;

	static TextField exptDirTF;

	static Color bgCol=new Color(255, 255, 240);

	static Choice runChoice, graphVarChoice, paramsChoice, mapVarChoice;

	JCheckBox annualCB, paramCB, globeCB;
	static JCheckBox changeCB;

	static ParamSetting ps;
	static Vector paramsToChange;
	static AdjstblParam ap;
	int paramIndex=0, oldParamIndex=0;

	Scrollbar depthS, paramS;

	TextField paramTF;

	public static Image mapIm, keyIm;

	public static Panel graphPanel, mapPanel
		, keyPanel=new Panel(){
				public void paint(Graphics g){
					super.paint(g);
					if(g!=null && keyIm!=null)g.drawImage(keyIm, 0, 0, null);
				}
			}
		, mapControlsPanel=new Panel()
		, ocean3DPanel=new Panel();

	public static Panel paramsPanel=new Panel();

  static Frame frame=new Frame(), chooserFrame;

	Frame paramFrame=new Frame();
	
	static TimeSeriesFrame tSF;

	int mapW0=360, mapW=mapW0, mapH0=360, mapH=mapH0
		, mapX=360
		, cellW=mapW/iMax, cellH=mapH/jMax
		, keyW=60, keyH=360;

	static TextDrawPanel graphTitleTDP=new TextDrawPanel(), depthTDP=new TextDrawPanel()
		, keyTDP=new TextDrawPanel(), mapTDP=new TextDrawPanel()
		, graphTDP=new TextDrawPanel(), graphTimeTDP=new TextDrawPanel()
		, paramTDP=new TextDrawPanel()
		, exptTDP=new TextDrawPanel()
		, info1TDP=new TextDrawPanel(), info2TDP=new TextDrawPanel()
		;

	Process process;

	static String os, exptTitle
		, info1S="<b>Experiment data folder selection</b><br/>"
        +"Select the folder containing the genie output data (either a folder containing genie.exe or a time-stamped results folder)."
				+"<br/><br/><b>Current folder:<br/>"
		, info2S="<b>Comparison data folder selection</b><br/>"
        +"Select the folder containing the genie output data that you wish to "
        +"compare with current results."
				+"<br/><b>Comparison folder:<br/>";

	public static boolean isWindows;

	public static int graphVarIdx=0, mapVarIdx=0, nKLayers=1, kLayer=0, nLLayers=1, lLayer=0
		, nOuterLayers=0, outerLayer=1, nInitialMaps=0, dirn=0, depth=7, sDepth=7
	;

	public static double scale=1., offset=0., mapValP, mapDPP;
	public static int mapValP10;

	public static boolean reading=false, isOcean=true
    , isOcean3D=true, isLand=false, isAll=false, paused=false
		, saveData=false, plotData=true
		, setChange=false;

	public static String[] runTags, graphVarTags, mapVarTags;
	public static String defnXMLFile, jobXMLFile
		, graphVarFile, mapVarFile
		, graphTitle="", graphUnits="", mapTitle="", mapUnits="";

	StringBuffer paramHTMLSB=new StringBuffer();

	Thread timer;
	int delay=1100;

	Graphics g;

	GraphPlotter graphPlotter;

	GlobePlotter globePlotter;

	HTMLFiles hf;

	static final int nMapColrs=200, nMapColrs1=nMapColrs-1;

	Color seaBedCol=new Color(128, 128, 128), greenCol=new Color(0, 110, 0)
		, blueCol=new Color(0, 0, 130), grey=new Color(220, 220, 220)
		, rainbowCols[]=new Color[nMapColrs]
		;

	public BasicStroke bs=new BasicStroke(3);

	public static GenieGUI genieGUI;

	public static MenuBar mb=new MenuBar();
	public static Menu fileMenu=new Menu("File"), helpMenu=new Menu("Help")
    , exptMenu=new Menu("Experiment"), runMenu=new Menu("Run")
    , recentExptsMenu=new Menu("Recent experiments");
  public static MenuItem pauseMI=new MenuItem("Pause")
    , paramsMI, inputMI
		, statusMI=new MenuItem("Status")
    , contMI=new MenuItem("Continue")
    , stopMI=new MenuItem("Stop");
	public static int menuRun;

  public static void main(String[] args) {

    genieGUI = new GenieGUI();

    genieGUI.isStandalone = true;
    frame.addWindowListener(new WindowAdapter() { 
			public void windowClosing(WindowEvent e) { System.exit(0); } 
		});
    frame.setTitle("Genie interface");
    frame.add(genieGUI, BorderLayout.CENTER);

		MenuItem helpMI=new MenuItem("Help", new MenuShortcut((int)'H'));
		helpMI.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
        try{
          String pageS="./help.html";
          if(isWindows){
            if (os.equals("Windows 95") || os.equals("Windows 98") || os.equals("Windows ME")) {
              Runtime.getRuntime().exec("start " + pageS);
            } else {
              Runtime.getRuntime().exec("cmd /c start \"name\" \"" + pageS + "\"");
            }
          }	else Runtime.getRuntime().exec("open "+pageS);
        }catch(Exception e1){}
			}
		});
		helpMenu.add(helpMI);
		
		MenuItem aboutMI=new MenuItem("About", new MenuShortcut((int)'A'));
		aboutMI.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Utils.showAbout("images/Logos.png", " ");
				Utils.aboutTDP.setLayout(null);
				int y0=30;
				LinkButton l1B=new LinkButton("Designed and built as part of the ENGAGE programme <b>http://www.engage.ac.uk", isWindows, os);
				l1B.setBackground(bgCol);
				l1B.setLocation(20,y0+4);
				Utils.aboutTDP.add(l1B, null);
				LinkButton l0B=new LinkButton("by the ALADDIN2 project team:", isWindows, os);
				l0B.setBackground(bgCol);
				l0B.setLocation(20,y0+24);
				Utils.aboutTDP.add(l0B, null);
				LinkButton l2B=new LinkButton("  Martin Johnson, Sudipta Goswami   <b>http://www.uea.ac.uk", isWindows, os);
				l2B.setBackground(bgCol);
				l2B.setLocation(20,y0+44);
				Utils.aboutTDP.add(l2B, null);
				LinkButton l3B=new LinkButton("  Will Rawes, Simon Mueller, Neil Edwards   <b>http://www.open.ac.uk", isWindows, os);
				l3B.setBackground(bgCol);
				l3B.setLocation(20,y0+64);
				Utils.aboutTDP.add(l3B, null);
				LinkButton l5B=new LinkButton("  Gethin Williams   <b>http://www.bristol.ac.uk", isWindows, os);
				l5B.setBackground(bgCol);
				l5B.setLocation(20,y0+84);
				Utils.aboutTDP.add(l5B, null);
				LinkButton l4B=new LinkButton("  Andrew Price   <b>http://www.soton.ac.uk", isWindows, os);
				l4B.setBackground(bgCol);
				l4B.setLocation(20,y0+104);
				Utils.aboutTDP.add(l4B, null);
				LinkButton l6B=new LinkButton("Using the GENIE earth-system model   <b>http://www.genie.ac.uk", isWindows, os);
				l6B.setBackground(bgCol);
				l6B.setLocation(20,y0+124);
				Utils.aboutTDP.add(l6B, null);
				LinkButton l7B=new LinkButton("For further information see   <b>http://researchpages.net/genie/aladdin", isWindows, os);
				l7B.setBackground(bgCol);
				l7B.setLocation(20,y0+144);
				Utils.aboutTDP.add(l7B, null);
			}
		});
		helpMenu.add(aboutMI);

		mb.add(fileMenu);
		mb.add(exptMenu);
		mb.add(runMenu);
		mb.add(helpMenu);

		frame.setMenuBar(mb);

    frame.setSize(minWidth, minHeight);
    Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
    frame.setLocation((d.width - frame.getSize().width) / 2, (d.height - frame.getSize().height) / 2);
    frame.setVisible(true);

    genieGUI.init();

		genieGUI.repaint();

    frame.addComponentListener(new ComponentAdapter() { 
			public void componentResized(ComponentEvent e) {
				genieGUI.resize(frame.getWidth(), frame.getHeight());
			}
		});
  }

	public void paint(Graphics g){
		super.paint(g);
	}

	public void resize(int w, int h){
		if(w<minWidth)w=minWidth;
		if(h<minHeight)h=minHeight;
		mapH=mapH0+h-600;
		int d=365+mapH+keyW+5-w;
		if(d>0)mapH-=d;
		mapW=mapH;
		cellW=mapW/iMax;
		cellH=mapH/jMax;
		mapW=cellW*iMax;
		mapH=mapW;
		mapX=w-(mapW+keyW+5)-15;
		d=mapX-365;
		if(d>0)graphPanel.setBounds(10,230, 350+d,300);
		else graphPanel.setBounds(10,230, 350,300);
		graphPlotter=new GraphPlotter(graphPanel);
		mapIm=createImage(mapW, mapH);
		mapPanel.setBounds(mapX,100, mapW+keyW+5,mapH+25+80);
		mapControlsPanel.setLocation(20, mapH+35);
		frame.setSize(w, h);
    drawGraph();
		drawMap();
	}

  public void init(){
		os = System.getProperty("os.name");
		isWindows=(os.indexOf("Windows")>-1);

		setLayout(null);

		GenieGUI genieGUI = new GenieGUI();

		hf=new HTMLFiles();

		for(int i=0;i<nMapColrs;i++)rainbowCols[i]=Colours.setRainbowCol((double)i, 0., nMapColrs1);

		readInitialData();

    nCoastlines=0;
    for (int i = 0; i < iMax; i++) {
      for (int j = 0; j < jMax1; j++) {
        if(land[i][j]!=land[(i+1)%iMax][j]){
          clStartX[nCoastlines]=(i+1)*10;
          clStartY[nCoastlines]=j*10;
          clEndX[nCoastlines]=(i+1)*10;
          clEndY[nCoastlines]=(j+1)*10;
          nCoastlines++;
        }
        if(land[i][j]!=land[i][(j+1)%jMax]){
          clStartX[nCoastlines]=i*10;
          clStartY[nCoastlines]=(j+1)*10;
          clEndX[nCoastlines]=(i+1)*10;
          clEndY[nCoastlines]=(j+1)*10;
          nCoastlines++;
        }
      }
    }

		graphPanel=new Panel(){
			public void paint(Graphics g){
				super.paint(g);
				if(g!=null && graphPlotter.graphIm!=null)g.drawImage(graphPlotter.graphIm, 0, 0, null);
			}
		};
		graphPanel.setBounds(10,230, 350,298);

		graphPlotter=new GraphPlotter(graphPanel);

		setBackground(graphPlotter.bgCol);

		globePlotter=new GlobePlotter(bs, greenCol, nCoastlines
			, clStartX, clStartY, clEndX, clEndY);

//	  String XMLContent=GenieFileReader.getXMLContents("runData.xml");
	  String XMLContent="<run>"
			+"	<name>ENGAGEtest</name>	<outputDir>.</outputDir>"
			+"	<runDir>ENGAGEtest</runDir>"
			+"	<initialDir>ENGAGEtest</initialDir>	<definitionXMLFile>definition.xml</definitionXMLFile>"
			+"	<jobXMLFile>ENGAGEtest.xml</jobXMLFile>	<graphVarFile>graphVarData.xml</graphVarFile>"
			+"	<mapVarFile>mapVarData.xml</mapVarFile>	<RUNTIME_ROOT>../../genie</RUNTIME_ROOT>"
			+"	<RUNTIME_OUTDIR>.</RUNTIME_OUTDIR>"
			+"</run>"	/// Default runData file settings
		; 

		runTags=XMLReader.getTags(XMLContent, "run");

		MenuItem options=new MenuItem("Open experiment/data folder...", new MenuShortcut((int)'O'));
		options.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
        setChange=false;
				setUpRunData();
			}
		});
		fileMenu.add(options);

    setPrevRunMenu(Utils.readFile("previousRuns.dat"));
    
/*		Menu exptsMenu=new Menu("Experiments");
		for(int j=0;j<runTags.length;j++){
			menuRun=j;
			MenuItem options2=new MenuItem(XMLReader.getTag(runTags[j], "name"));//, new MenuShortcut((int)((char)('a'+j))));
			options2.addActionListener(new ActionListener() {
				int lMenuRun=menuRun;
				public void actionPerformed(ActionEvent e) {
					readRunDataXML(lMenuRun);
					setParamPanel();
				}
			});
			exptsMenu.add(options2);
		}
		fileMenu.add(exptsMenu);//*/
		
		fileMenu.add(recentExptsMenu);

		MenuItem exit=new MenuItem("Exit", new MenuShortcut((int)'X'));
		exit.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {System.exit(0);}
		});
		fileMenu.add(exit);
		fileMenu.insertSeparator(2);

		paramsMI=new MenuItem("Adjust parameters", new MenuShortcut((int)'P'));
		paramsMI.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				setParamControl();
				paramFrame.setVisible(true);
			}
		});
		exptMenu.add(paramsMI);

		inputMI=new MenuItem("Adjust input data", new MenuShortcut((int)'I'));
		inputMI.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
        tSF.genieRunDir=genieRunDir;
				tSF.setVisible(true);
			}
		});
		exptMenu.add(inputMI);


///////////////
		statusMI.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
        Utils.showFile(genieRunDir+"/genie.out");
			}
    });
		runMenu.add(statusMI);

		pauseMI.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				paused=true;
        Utils.writeFile(genieRunDir+"/genie_pause", " ");
        pauseMI.setEnabled(false);
        contMI.setEnabled(true);
			}
    });
    pauseMI.setEnabled(false);
		runMenu.add(pauseMI);

		contMI.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				paused=false;
        Utils.writeFile(genieRunDir+"/genie_continue", " ");
        pauseMI.setEnabled(true);
        contMI.setEnabled(false);
			}
    });
    contMI.setEnabled(false);
		runMenu.add(contMI);

    stopMI.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
        try{
          if(isWindows){
            if(process!=null){
              process.destroy();
System.out.println("Attempt to stop: "+process);
            }
            Utils.showMessage("<br/><br/>Please close the appropriate DOS cmd window if necessary to stop the run.");
            Utils.aboutOKB.grabFocus();
          }else{
            File f=new File(genieRunDir+"/genie_active");
            f.delete();
          }
        }catch(Exception e2){
          System.out.println("Minor error stopping run: "+e2);
        }
				saveData=true;
				stopMI.setEnabled(false);
				pauseMI.setEnabled(false);
				contMI.setEnabled(false);
				startB.setEnabled(true);
			}
    });
		stopMI.setEnabled(false);
		runMenu.add(stopMI);

///////////

		ps=new ParamSetting();

		tSF=new TimeSeriesFrame(graphPlotter.bgCol);

    paramFrame.addWindowListener(new WindowAdapter() { 
			public void windowClosing(WindowEvent e) { paramFrame.setVisible(false); }
		});
    paramFrame.setTitle("Genie interface - Parameter adjustment");
    paramFrame.setSize(500, 380);
//    paramFrame.setLocation(30, 30);
    Utils.centreFrame(paramFrame);
    paramFrame.setVisible(false);

		paramsPanel=new Panel(){
			public void paint(Graphics g){
				super.paint(g);
				Graphics2D g2=(Graphics2D)g;
				g2.setColor(Color.black);
				g2.drawString("Parameter settings", 10,20);
				g2.drawString("Parameter settings", 11,20);
			}
		};
		paramsPanel.setLayout(null);
		paramsPanel.setVisible(true);//false);
		paramsPanel.setBackground(graphPlotter.bgCol);
    paramFrame.add(paramsPanel, BorderLayout.CENTER);

		paramsChoice=new Choice();
		paramsChoice.setBounds(10,30, 430,25);
		paramsChoice.addItemListener(new ItemListener() {
			public void itemStateChanged (ItemEvent e){
				setParamControl();
		}});
		paramsPanel.add(paramsChoice);

		readRunDataXML(0);
		setUpRunData();

		paramS=new Scrollbar(Scrollbar.HORIZONTAL, 7, 2, 0, 7+2);
		paramS.setBounds(10,65, 300,15);
		paramS.addAdjustmentListener(new AdjustmentListener(){
			public void adjustmentValueChanged(AdjustmentEvent e){
				ap=(AdjstblParam)(paramsToChange.get(paramIndex));
				ap.temp=""+(paramS.getValue()*ap.interval+ap.min);
				paramTF.setText(""+ap.temp);
		}});
		paramsPanel.add(paramS, null);

		paramCB=new JCheckBox("State");
		paramCB.setBounds(10,65, 300,25);
		paramCB.setFocusable(true);
		paramCB.setBackground(graphPlotter.bgCol);
		paramCB.setSelected(false);
		paramCB.addItemListener(new ItemListener() {
		public void itemStateChanged(ItemEvent e) {
			ap=(AdjstblParam)(paramsToChange.get(paramIndex));
			ap.temp=(paramCB.isSelected()?".true.":".false.");
		}});
		paramsPanel.add(paramCB, null);

		paramTF=new TextField();
		paramTF.setBounds(320,60, 120,25);
		paramTF.addKeyListener(new KeyListener(){
			public void keyReleased(KeyEvent e){
				ap=(AdjstblParam)(paramsToChange.get(paramIndex));
/*				double v;   // This could could be used to check for valid number input but doesn't allow text..
				try{					// ... so is commented out
					v=Double.valueOf(paramTF.getText());
				}catch(Exception e2){
					paramTDP.setText("Please input a valid number.");
					return;
				}//*/
				ap.temp=(paramTF.getText());
				paramTDP.setText("<b>Units</b>: "+ap.units+"<br/><b>Description</b>: "+ap.description+"<br/><b>Namelist</b>: "+ap.namelist);
			}
			public void keyTyped(KeyEvent e){}
			public void keyPressed(KeyEvent e){
		}});
		paramsPanel.add(paramTF);

		paramTDP.setBounds(10,90, 470,200);
		paramTDP.setExtraTagsEnabled(false);
		paramTDP.setBackground(graphPlotter.bgCol);
		paramTDP.setFont("SansSerif", Font.PLAIN, 12);
		paramsPanel.add(paramTDP, null);

		paramsDefB=new KButton("Set default values", this, true); 
    paramsDefB.setLocation(20,310);
    paramsDefB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				paramsDefB.down=true;
				paramsDefB.paint(paramsDefB.getGraphics());
				setParamPanel();
				setParamControl();
				paramsDefB.down=false;
				paramsDefB.paint(paramsDefB.getGraphics());
			}
    });
		paramsPanel.add(paramsDefB, null);

		paramsCanB=new KButton("Cancel", this, true); 
    paramsCanB.setLocation(220,310);
    paramsCanB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				plotData=true;
				paramsCanB.down=true;
				paramsCanB.paint(paramsCanB.getGraphics());
				paramFrame.setVisible(false);
				paramsCanB.down=false;
				paramsCanB.paint(paramsCanB.getGraphics());
			}
    });
		paramsPanel.add(paramsCanB, null);

		paramsOKB=new KButton("OK", this, true); 
    paramsOKB.setLocation(340,310);
    paramsOKB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				plotData=true;
				paramsOKB.down=true;
				paramsOKB.paint(paramsOKB.getGraphics());
				paramHTMLSB=new StringBuffer();
				paramHTMLSB.append("<html>\n<head>\n<link href='../style.css' rel='stylesheet' type='text/css'/>\n"
					+"</head>\n<body>\n"//<<HTMLlinks<<"\n"
					+"<h5>Changed parameter values:"
					+"</h5>"
					+"<table cellspacing='0' id='tbl'>\n"
					+"<tr><td><b>Name</b></td><td><b>Description</b></td><td><b>Value</b></td></tr>\n"
				);
				setParamVals();
//				Utils.writeFile(folder+"/"+filename, paramHTMLSB.toString());
				ps.runXSLT(defnXMLFile, jobXMLFile, genieRunDir, paramsToChange.toArray());
				paramFrame.setVisible(false);
				paramsOKB.down=false;
				paramsOKB.paint(paramsOKB.getGraphics());
			}
    });
		paramsPanel.add(paramsOKB, null);

		add(graphPanel, null);

		mapPanel=new Panel(){
			public void paint(Graphics g){
				super.paint(g);
				if(g!=null && mapIm!=null)g.drawImage(mapIm, keyW+2,26, null);
			}
		};
		mapPanel.setBounds(mapX,100, mapW+keyW+5,mapH+25+80);
		mapPanel.addMouseListener(new MapML());
		mapPanel.addMouseMotionListener(new MapMML());
		mapPanel.setLayout(null);
		add(mapPanel, null);

		setParamPanel();

System.out.println("run tags content:");
		for(int i=0;i<runTags.length;i++)System.out.println(runTags[i]);

		exptTDP.setBounds(10,3, 775,19);
		exptTDP.setExtraTagsEnabled(false);
		exptTDP.setBackground(graphPlotter.bgCol);
		exptTDP.setFont("SansSerif", Font.PLAIN, 12);
		add(exptTDP, null);
/*		runChoice=new Choice();
		runChoice.setFont(new Font("SanSerif", Font.PLAIN, 11));
		runChoice.setBounds(10,10, 250,25);
		for(int j=0;j<runTags.length;j++)runChoice.add(XMLReader.getTag(runTags[j], "name"));
		runChoice.addItemListener(new ItemListener() {
			public void itemStateChanged (ItemEvent e){
				readRunDataXML(runChoice.getSelectedIndex());
				setParamPanel();
System.out.println("genie Output & command:"+genieRunDir+" "+genieCommand);
		}});
		add(runChoice, null);//*/

/*		paramsB=new KButton("Adjust parameters", this, true); 
    paramsB.setLocation(10,55);
    paramsB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				plotData=false;
				paramsB.down=true;
				paramsB.paint(paramsB.getGraphics());
				setParamControl();
				paramFrame.setVisible(true);
				paramsB.down=false;
				paramsB.paint(paramsB.getGraphics());
			}
    });
		add(paramsB, null);

		emsB=new KButton("Adjust input file", this, true); 
    emsB.setLocation(180,55);
    emsB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				emsB.down=true;
				emsB.repaint();

				tSF.setVisible(true);

				emsB.down=false;
				emsB.repaint();
			}
    });
		add(emsB, null);//*/

		startB=new KButton("Start", this, true); 
    startB.setLocation(10,40);
		startB.setFocusable(true);
    startB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				plotData=true;
				startB.down=true;
				startB.repaint();

				try{
					String envp[]=null;
System.out.println("os:"+os);
					if(isWindows){
System.out.print("Running:genie.bat in "+genieRunDir);
						Utils.writeFile(genieRunDir+"/genie.bat", "genie.exe>genie.out\nexit\n");
						if (os.equals("Windows 95") || os.equals("Windows 98") || os.equals("Windows ME")) {
							process=Runtime.getRuntime().exec("start " + "genie.bat", envp, new File(genieRunDir));
						} else {
							process=Runtime.getRuntime().exec("cmd /c start \"name\" \"" + "genie.bat" + "\"", envp, new File(genieRunDir) );
							stopMI.setEnabled(true);
							pauseMI.setEnabled(false);
						}
            stopMI.setEnabled(true);
            pauseMI.setEnabled(false);
            contMI.setEnabled(false);
					}	else {//if (os.indexOf("Mac")>=0 || os.indexOf("inux")>=0){
						//Utils.copyFile("./genie_wrapper.sh", genieRunDir+"/genie_wrapper.sh");
System.out.println("Running: genie_wrapper.sh in "+genieRunDir);
						process=Runtime.getRuntime().exec("./genie_wrapper.sh", envp, new File(genieRunDir));
            stopMI.setEnabled(true);
            pauseMI.setEnabled(true);
            contMI.setEnabled(true);
					} 
					startB.setEnabled(false);
				} catch (Exception e2){		
					System.out.println("Error running genie.exe: "+e2);
				}

				hf.createHTMLFiles(paramHTMLSB.toString(), genieRunDir);
				saveTestXMLFile();

				startB.down=false;
				startB.repaint();
			}
    });
		add(startB, null);

		graphTDP.setBounds(10,85, graphPanel.getWidth(),36);
		graphTDP.setExtraTagsEnabled(false);
		graphTDP.setBackground(graphPlotter.bgCol);
		graphTDP.setFont("SansSerif", Font.PLAIN, 12);
		add(graphTDP, null);

/*		pauseB=new KButton("Pause", this, true); 
    pauseB.setLocation(10,100);
    pauseB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				pauseB.down=true;
				pauseB.paint(pauseB.getGraphics());
				paused=!paused;
				if(paused){
					Utils.writeFile(genieRunDir+"/genie_pause", " ");
					pauseB.setText("Continue");
//					pauseB.upIm=new ImageIcon("images/continueUp.png").getImage();
				}else{
					Utils.writeFile(genieRunDir+"/genie_continue", " ");
					pauseB.setText("Pause");
//					pauseB.upIm=new ImageIcon("images/pauseUp.png").getImage();
				}
				pauseB.setSize(pauseB.upIm.getWidth(null),pauseB.upIm.getHeight(null));
				pauseB.down=false;
				pauseB.paint(pauseB.getGraphics());
			}
    });
		pauseB.setVisible(false);
		add(pauseB, null);

		stopB=new KButton("Stop", this, true); 
    stopB.setLocation(220,100);
    stopB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				stopB.down=true;
				stopB.paint(stopB.getGraphics());
        try{
          if(isWindows){
            if(process!=null){
              process.destroy();
System.out.println("Attempt to stop: "+process);
            }
            Utils.showMessage("<br/><br/>Please close the appropriate DOS cmd window if necessary to stop the run.");
            Utils.aboutOKB.grabFocus();
          }else{
            File f=new File(genieRunDir+"/genie_active");
            f.delete();
          }
        }catch(Exception e2){
          System.out.println("Minor error stopping run: "+e2);
        }
				saveData=true;
				stopB.setVisible(false);
				pauseB.setVisible(false);
				startB.setVisible(true);
				stopB.down=false;
				stopB.paint(stopB.getGraphics());
			}
    });
		stopB.setVisible(false);
		add(stopB, null);//*/

/// Read graph variable xml file
		graphVarChoice=new Choice();
		graphVarChoice.setFont(new Font("SanSerif", Font.PLAIN, 11));
		graphVarChoice.setBounds(10,140, 250,25);
		readGraphVarDataXML(graphVarFile);
		graphVarChoice.addItemListener(new ItemListener() {
			public void itemStateChanged (ItemEvent e){
				graphVarIdx=graphVarChoice.getSelectedIndex();
				String tagText=graphVarTags[graphVarIdx];
				timeSFile=Utils.removeWhitespace(XMLReader.getTag(tagText, "file"));
				colNum=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "col"))-1;
				nCols=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "nCols"));
				nHeaderLines=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "nHeaderLines"));
				graphTitle=XMLReader.getTag(tagText, "name");
				graphUnits=XMLReader.getTag(tagText, "units");
				readTimeSeriesFile(genieRunDir+File.separator+timeSFile);
System.out.println("ts data:"+timeSFile+" "+colNum+" "+nCols+" "+nHeaderLines);
		}});
		graphVarIdx=0;
		add(graphVarChoice, null);

		annualCB=new JCheckBox("Annual averages");
		annualCB.setBounds(10,160, 250,25);
		annualCB.setBackground(graphPlotter.bgCol);
		annualCB.setFocusable(true);
		annualCB.setSelected(false);
		annualCB.addItemListener(new ItemListener() {
		public void itemStateChanged(ItemEvent e) {
	    drawGraph();
		}});
		add(annualCB, null);
		annualCB.repaint();

		graphTitleTDP.setBounds(40,200, 320,25);
		graphTitleTDP.setExtraTagsEnabled(false);
		graphTitleTDP.setBackground(graphPlotter.bgCol);
		graphTitleTDP.setFont("SansSerif", Font.PLAIN, 12);
		graphTitleTDP.setText(graphTitle);//[graphVarIdx]);
		add(graphTitleTDP, null);

		graphTimeTDP.setBounds(graphPanel.getLocation().x+180
			,graphPanel.getLocation().y+graphPanel.getHeight()
			, 120,25);
		graphTimeTDP.setExtraTagsEnabled(false);
		graphTimeTDP.setBackground(graphPlotter.bgCol);
		graphTimeTDP.setFont("SansSerif", Font.PLAIN, 12);
		graphTimeTDP.setText(graphTitle);//[graphVarIdx]);
		add(graphTimeTDP, null);
		graphTimeTDP.setText("time / years");

		mapTDP.setBounds(mapPanel.getLocation().x,mapPanel.getLocation().y-50
      , mapW, 50);
		mapTDP.setExtraTagsEnabled(false);
		mapTDP.setBackground(graphPlotter.bgCol);
		mapTDP.setFont("SansSerif", Font.PLAIN, 12);
		add(mapTDP, null);

		keyPanel.setBounds(0, 25, keyW, keyH);
		mapPanel.add(keyPanel, null);

		mapControlsPanel.setBounds(20,mapH+35, mapW0+50,100);
		mapControlsPanel.setLayout(null);
		mapPanel.add(mapControlsPanel, null);

		keyTDP.setBackground(graphPlotter.bgCol);
		keyTDP.setFont("SansSerif", Font.PLAIN, 12);
		mapPanel.add(keyTDP, null);

		keyTDP.setBounds(0,0, 300,19);
		keyTDP.setExtraTagsEnabled(false);
		keyTDP.setBackground(graphPlotter.bgCol);
		keyTDP.setFont("SansSerif", Font.PLAIN, 12);
		mapPanel.add(keyTDP, null);

/// Read map variable xml file
		mapVarChoice=new Choice();
		mapVarChoice.setFont(new Font("SanSerif", Font.PLAIN, 11));
		mapVarChoice.setBounds(5,0, 210,30);
		readMapVarDataXML(mapVarFile);
		mapVarChoice.addItemListener(new ItemListener() {
			public void itemStateChanged (ItemEvent e){
				mapVarIdx=mapVarChoice.getSelectedIndex();
				setMapFile(mapVarTags[mapVarIdx]);
				readInitialMapVal(mapVarIdx);
				drawMap(readMapFile(genieRunDir+File.separator+mapFile, changeCB.isSelected()));
		}});
		mapVarIdx=0;
		mapControlsPanel.add(mapVarChoice, null);

		ocean3DPanel.setLayout(null);
		ocean3DPanel.setBounds(5,30, 220,50);
		mapControlsPanel.add(ocean3DPanel, null);

		depthS=new Scrollbar(Scrollbar.HORIZONTAL, 7, 2, 0, 7+2);
		depthS.setBounds(72,0, 140,15);
		depthS.addAdjustmentListener(new AdjustmentListener(){
			public void adjustmentValueChanged(AdjustmentEvent e){
				int oldDepth=sDepth;
				sDepth=depthS.getValue();
				if(oldDepth==sDepth)return;
				depthTDP.setText("Layer: "+sDepth);
				setMapFile(mapVarTags[mapVarIdx]);
				readInitialMapVal(mapVarIdx);
				drawMap(readMapFile(genieRunDir+File.separator+mapFile, changeCB.isSelected()));
		}});
		depthS.addKeyListener(new KeyAdapter(){
			public void keyPressed(KeyEvent e){
				drawMap(readMapFile(genieRunDir+File.separator+mapFile, changeCB.isSelected()));
			}
		});
		depthS.addMouseListener(new MouseAdapter(){
			public void mouseReleased(MouseEvent e){
				drawMap(readMapFile(genieRunDir+File.separator+mapFile, changeCB.isSelected()));
			}
		});
		ocean3DPanel.add(depthS, null);

		depthTDP.setBounds(0,0, 70, 20);
		depthTDP.setExtraTagsEnabled(false);
		depthTDP.setBackground(graphPlotter.bgCol);
		depthTDP.setFont("SansSerif", Font.PLAIN, 12);
		depthTDP.setText("Layer: "+sDepth);
		ocean3DPanel.add(depthTDP, null);

		changeB=new KButton("Compare", this, true); 
		changeB.setBounds(240,0, 140,25);
		changeB.setFocusable(true);
    changeB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				changeB.down=true;
				changeB.repaint();
				
        setChange=true;
        setUpRunData();
        
				changeB.down=false;
				changeB.repaint();
			}
    });
		mapControlsPanel.add(changeB, null);
		changeB.repaint();

		globeCB=new JCheckBox("Show globe");
		globeCB.setFocusable(true);
		globeCB.setBounds(240,25, 240,25);
		globeCB.setBackground(graphPlotter.bgCol);
		globeCB.setSelected(false);
		globeCB.addItemListener(new ItemListener() {
		public void itemStateChanged(ItemEvent e) {
			drawMap();
		}});
		mapControlsPanel.add(globeCB, null);
		globeCB.repaint();

		timer=new Thread(this);
		timer.setPriority(Thread.MIN_PRIORITY+2);
		timer.start();
	}

	int c=0;
	public void run(){
		while(timer!=null){
			try{
				timer.sleep( 1700 );
			}catch(Exception e){}
			if(plotData){
		    readTimeSeriesFile(genieRunDir+File.separator+timeSFile);
				drawMap(readMapFile(genieRunDir+File.separator+mapFile, changeCB.isSelected()));
			}
			if(saveData){
				saveData=false;
				saveGraphData();
				saveMapData();
			}
		}
	}

	public static void setParamPanel(){
		try{
			paramsToChange=ps.getParams(defnXMLFile, jobXMLFile, genieRunDir);
System.out.println("Defn, job, rundir, paramsToChange size:\n"+defnXMLFile+"\n"+ jobXMLFile+"\n"+ genieRunDir+"\n"+paramsToChange.size());
			paramsChoice.removeAll();
			for(int i=0;i<paramsToChange.size();i++){
				ap=(AdjstblParam)(paramsToChange.get(i));
				paramsChoice.add( ap.name );
			}

			if(paramsToChange.size()==0)Utils.showMessage("There are no exposed parameters for this experiment.");
			else paramsChoice.select(0);
		}catch(Exception e){
			Utils.showMessage("Error getting parameter data from<br/>Definition.xml file:"+defnXMLFile
				+",<br/>Experiment file: "+jobXMLFile+"<br/>for folder:<br/>"+genieRunDir);
		}
	}

	public void setParamControl(){
		try{
			paramIndex=paramsChoice.getSelectedIndex();
			ap=(AdjstblParam)(paramsToChange.get(oldParamIndex));
			if(ap.temp!=null)ap.value=ap.temp;
			oldParamIndex=paramIndex;
			ap=(AdjstblParam)(paramsToChange.get(paramIndex));
			ap.temp=ap.value;
			if(ap.temp==null){
System.out.println("\n!!! null value for "+ap+" "+paramIndex);
return;
			}
			paramTDP.setText("<b>Units</b>: "+ap.units+"<br/><b>Description</b>: "+ap.description+"<br/><b>Namelist</b>: "+ap.namelist);
			paramCB.setVisible(false);
			paramS.setVisible(false);
			paramTF.setVisible(false);
			if(ap.temp.equalsIgnoreCase(".true.") || ap.temp.equalsIgnoreCase(".false.")){
				paramCB.setSelected(ap.temp.equalsIgnoreCase(".true."));
				paramCB.setVisible(true);
			}else if(ap.interval==0){
				paramTF.setVisible(true);
			}else{
				paramS.setMinimum(0);//(int)(ap.min/ap.interval));
				paramS.setMaximum((int)((ap.max-ap.min)/ap.interval)+2);
				paramS.setValue((int)( (Double.valueOf(ap.temp)-ap.min)/ap.interval ));
				paramS.setVisible(true);
				paramTF.setVisible(true);
			}
			paramTF.setText(""+ap.temp);
		}catch(Exception e){
			Utils.showMessage("Error setting adjustable parameters.<br/>"
				);//+paramsToChange.size());
		}
	}

	public void setParamVals(){
		for(int i=0;i<paramsToChange.size();i++){
			ap=(AdjstblParam)(paramsToChange.get(i));
			if(ap.temp!=null){
				ap.value=ap.temp;
				paramHTMLSB.append("<tr><td>"+ap.name+"</td><td>"+ap.description+"</td><td>"+ap.value+" "+ap.units+"</td></tr>\n");
			}
		}
	}

	public void saveTestXMLFile(){
		Utils.copyFile("nmlBuilding/jobTMP.xml", hf.runFolder+File.separator+"test.xml");
	}

	public void saveGraphData(){
		for(int i=0;i<graphVarTags.length;i++){
			String graphFileName=Utils.removeWhitespace(XMLReader.getTag(graphVarTags[i], "file"));
			Utils.copyFile(genieRunDir+File.separator+graphFileName
				, hf.runFolder+File.separator+graphFileName.substring(0,graphFileName.indexOf("/"))
				, hf.runFolder+File.separator+graphFileName);
		}
	}
	
	public void saveMapData(){
		for(int i=0;i<mapVarTags.length;i++){
			String mapFileName=Utils.removeWhitespace(XMLReader.getTag(mapVarTags[i], "file"));
			Utils.copyFile(genieRunDir+File.separator+mapFileName
				, hf.runFolder+File.separator+mapFileName.substring(0,mapFileName.indexOf("/")) 
				, hf.runFolder+File.separator+mapFileName);
		}
	}
	
/////////////// Draws graph of timeVal[i] vs. selected variable
	public void drawGraph(){
	  if(annualCB.isSelected()){
			int nPs=0;
			int tv0=(int)timeVal[0], tv=tv0;
			while(nPs<nReadings && tv==tv0){
				nPs++;
				tv=(int)timeVal[nPs];
			}
			nPs++;
System.out.println("no. points in average:"+nPs);
	    xScale=nPs;
	    xMin=1;
	    nYVals=(int)(nReadings/(double)nPs);
	    if(true){
		    yMin=1e200;
		    yMax=-1e200;
        for(int i=0;i<nYVals;i++){
          int j=i*nPs;
					yVals[i]=0;
					for(int k=0;k<nPs;k++)yVals[i]+=graphVal[j+k];
					yVals[i]/=nPs;
          if(yVals[i]>yMax)yMax=yVals[i];
          if(yVals[i]<yMin)yMin=yVals[i];
        }
	    }
	    xMax=nYVals;
	  }else{
	    xScale=1;
	    yMin=1e200;
	    yMax=-1e200;
	    xMin=timeVal[0];
	    if(nReadings>0)xMax=timeVal[nReadings-1];
	    nYVals=nReadings;
      for(int i=0;i<nYVals;i++){
        yVals[i]=graphVal[i];
        if(yVals[i]>yMax)yMax=yVals[i];
        if(yVals[i]<yMin)yMin=yVals[i];
      }
	  }
	  if(xMax==xMin)xMax=xMin+1;
	
	  if(yMin==yMax){
	    yMin-=1;
	    yMax+=1;
	  }
	
		if(nYVals<2)graphPlotter.plotBlank("Waiting for data from "+genieRunDir+File.separator+timeSFile);
		else graphPlotter.plotGraph(graphPanel.getGraphics(), timeVal, xScale, xMin, xMax, yVals, nYVals, yMin, yMax);
		graphTitleTDP.setText(graphTitle+" / "+graphPlotter.graphTitlePowerS+" "+graphUnits);
	}
//////////////////////

	void setRunVariables(){
		genieRunDir=exptDirTF.getText();
	}

  static void setPrevRunMenu(String pRS){
    int op=0, p=pRS.indexOf("\n"), i=0;
    recentExptsMenu.removeAll();
    while(p>=0 && op>=0){
      final MenuItem runMI=new MenuItem(pRS.substring(op, p));
      runMI.addActionListener(new ActionListener() {
        MenuItem mi=runMI;
        public void actionPerformed(ActionEvent e) {
System.out.println(mi.getLabel());
          genieRunDir=mi.getLabel();
          getXMLFilenames();
          showBtns();
          if(graphVarChoice!=null && mapVarChoice!=null){
            readGraphVarDataXML(graphVarFile);
            readMapVarDataXML(mapVarFile);
          }
        }
      });
      recentExptsMenu.add(runMI);
      op=p+1;
      p=pRS.indexOf("\n", p+1);
    }
  }

	static void setUpRunData(){
		plotData=false;
		if(chooserFrame==null){
			chooserFrame=new Frame();
	    chooserFrame.addWindowListener(new WindowAdapter() { 
				public void windowClosing(WindowEvent e) { chooserFrame.setVisible(false); } 
			});
			chooserFrame.setLayout(null);
			chooserFrame.setResizable(false);
			chooserFrame.setBackground(bgCol);
	
	    chooserFrame.setSize(400,550);
	    Utils.centreFrame(chooserFrame);

			exptDirTF=new TextField();
/*			exptDirTF.setBounds(20,75, 360,27);
			chooserFrame.add(exptDirTF, null);//*/

			info1TDP.setBounds(20,40, 360, 130);
			info1TDP.setBackground(bgCol);
			chooserFrame.add(info1TDP, null);

			exptDirB=new KButton("Browse", null, true);
	    exptDirB.setLocation(200-40,190);
	    exptDirB.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent e) {
					exptDirB.down=true;
					exptDirB.repaint();

			    JFileChooser chooser = new JFileChooser();
					chooser.setCurrentDirectory(new File(genieRunDir));
					chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			    int returnVal = chooser.showOpenDialog(chooserFrame);
			    if(returnVal == JFileChooser.APPROVE_OPTION) {
						if(setChange){
              info1TDP.setText(info2S+chooser.getSelectedFile().getPath());
              info2TDP.setText("<b>"+chooser.getSelectedFile().getPath()
                +" (time = "+getCompareTime(chooser.getSelectedFile().getPath())
                +" years)"
                +"<br/><br/>and<br/><br/>"+genieRunDir+"</b>");
						}else info1TDP.setText(info1S+chooser.getSelectedFile().getPath());
						exptDirTF.setText(chooser.getSelectedFile().getPath());
			    }

					exptDirB.down=false;
					exptDirB.repaint();
				}
	    });
			chooserFrame.add(exptDirB, null);
			
      changeCB=new JCheckBox("Compare changes between data from:");
      changeCB.setBounds(20,250, 280,25);
      changeCB.setFocusable(true);
      changeCB.setBackground(bgCol);
      changeCB.setSelected(false);
			chooserFrame.add(changeCB,null);

			info2TDP.setBounds(20,280, 360,210);
			info2TDP.setBackground(bgCol);
			chooserFrame.add(info2TDP, null);

			cCanB=new KButton("Cancel", null, true);
	    cCanB.setLocation(100,500);
	    cCanB.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent e) {
					plotData=true;
					cCanB.down=true;
					cCanB.repaint();
					chooserFrame.setVisible(false);
					cCanB.down=false;
					cCanB.repaint();
				}
	    });
			chooserFrame.add(cCanB, null);

			cOKB=new KButton("OK", null, true);
	    cOKB.setLocation(250,500);
	    cOKB.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent e) {
					plotData=true;
					cOKB.down=true;
					cOKB.repaint();

          if(setChange){
            genieInitialDir=exptDirTF.getText();
          }else{
            genieRunDir=exptDirTF.getText();
            String s=Utils.readFile("previousRuns.dat");
            if(s.indexOf(genieRunDir)>=0)s=s.replace(genieRunDir+"\n", "");
            s=genieRunDir+"\n"+s;
            setPrevRunMenu(s);
            Utils.writeFile("previousRuns.dat", s);
          }
					
					getXMLFilenames();
					showBtns();
					chooserFrame.setVisible(false);
					if(graphVarChoice!=null && mapVarChoice!=null){
						readGraphVarDataXML(graphVarFile);
						readMapVarDataXML(mapVarFile);
					}

					cOKB.down=false;
					cOKB.repaint();
				}
	    });
			chooserFrame.add(cOKB, null);
		}
		if(setChange){
      exptDirTF.setText(genieInitialDir);
	    chooserFrame.setTitle("Comparison setup");
			info1TDP.setText(info2S+genieInitialDir);
			changeCB.setVisible(true);
			info2TDP.setText("<b>"+genieInitialDir
        +" (time = "+compareTime+" years)"
        +"<br/>and<br/>"+genieRunDir+"</b>");
			info2TDP.setVisible(true);
		}else{
      exptDirTF.setText(genieRunDir);
	    chooserFrame.setTitle("Experiment data folder selection");
			info1TDP.setText(info1S+genieRunDir);
			changeCB.setVisible(false);
			info2TDP.setVisible(false);
		}
    chooserFrame.setVisible(true);
	}

	public static void getXMLFilenames(){
System.out.println("Get XML filenames");
		File f1=new File(genieRunDir+"/archive");
		String[] l=f1.list();
		defnXMLFile="./definition.xml";
		jobXMLFile="./ENGAGEtest.xml";
		String runTimeRoot="../../genie";
		if(l!=null){
System.out.println(genieRunDir+"/archive contents:");
			for(int i=0;i<l.length;i++){
System.out.println(l[i]);
				if(l[i].indexOf(".xml")>0){
					if(l[i].equals("definition.xml"))defnXMLFile=genieRunDir+"/archive/"+l[i];
					else if(l[i].equals("runTimeRoot.xml")){
						String XMLContent=GenieFileReader.getXMLContents(genieRunDir+"/archive/runTimeRoot.xml");
						runTimeRoot=XMLReader.getTag(XMLContent, "RUNTIME_ROOT");
					}else jobXMLFile=genieRunDir+"/archive/"+l[i];
				}
			}
		}
System.out.println("RunTimeRoot:"+runTimeRoot);
		ps.runTimeRoot=runTimeRoot;
		File f2=new File(genieRunDir+"/aladdin");
		String[] l2=f2.list();
System.out.println(genieRunDir+"/aladdin contents:"+l2);
		graphVarFile="./graphVarData.xml";
		mapVarFile="./mapVarData.xml";
		if(l2!=null){
			for(int i=0;i<l2.length;i++){
System.out.println(l2[i]);
				if(l2[i].indexOf(".xml")>0){
					if(l2[i].equals("mapVarData.xml"))mapVarFile=genieRunDir+"/aladdin/"+l2[i];
					if(l2[i].equals("graphVarData.xml"))graphVarFile=genieRunDir+"/aladdin/"+l2[i];
				}
			}
		}
	}

	static void showBtns(){
		File f=new File(genieRunDir+"/genie.exe");
		int p=genieRunDir.lastIndexOf("/");
		exptTitle=genieRunDir.substring(p+1);
		tSF.setVisible(false);
		if(!f.exists()){
//			exptTDP.setText("<b>Data folder: "+exptTitle+"</b>");
		}else{
			setParamPanel();
//			exptTDP.setText("<b>Experiment: "+exptTitle+"</b>");
		}
		exptTDP.setText("<b>Data folder: </b>"+Utils.shortenName(genieRunDir, 100));
/*		if(paramsB==null)return;
		paramsB.setEnabled(f.exists());
		emsB.setEnabled(f.exists());//*/
		if(startB==null)return;
		paramsMI.setEnabled(f.exists());
		inputMI.setEnabled(f.exists());
		startB.setEnabled(f.exists());
//		pauseB.setVisible(false);
//		stopB.setVisible(false);
		pauseMI.setEnabled(false);
		contMI.setEnabled(false);
		stopMI.setEnabled(false);
	}

	void readRunDataXML(int runTagIdx){
			String tagText=runTags[runTagIdx];
			genieOutputDir=XMLReader.getTag(tagText, "outputDir");
			if(isWindows)genieOutputDir=genieOutputDir.replace("/cygdrive/c","c:");

			genieRunDir=genieOutputDir+"/"+XMLReader.getTag(tagText, "runDir");

			genieInitialDir=genieOutputDir+"/"+XMLReader.getTag(tagText, "initialDir");
			getXMLFilenames();// Get XML filenames from genieRunDir/archive folder
System.out.println("Graph & map var xml files:"+graphVarFile+"\n"+mapVarFile);
	
			ps.runTimeRoot=XMLReader.getTag(tagText, "RUNTIME_ROOT");
			ps.runTimeOutDir=XMLReader.getTag(tagText, "RUNTIME_OUTDIR");
	
			showBtns();

			if(graphVarChoice!=null && mapVarChoice!=null){
				readGraphVarDataXML(graphVarFile);
				readMapVarDataXML(mapVarFile);
			}
	}

	public static void readGraphVarDataXML(String varFile){
/// Read graph variable xml file
		try{
		  String XMLContent=GenieFileReader.getXMLContents(varFile);
			graphVarTags=XMLReader.getTags(XMLContent+"<end></end>", "variable"); //! end tag needed to fix tag reading error
System.out.println("graphVar tags content:");
			for(int i=0;i<graphVarTags.length;i++)System.out.println(graphVarTags[i]);
			graphVarIdx=0;
			String tagText=graphVarTags[0];
			timeSFile=Utils.removeWhitespace(XMLReader.getTag(tagText, "file"));
			colNum=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "col"))-1;
			nCols=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "nCols"));
			nHeaderLines=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "nHeaderLines"));
			graphTitle=XMLReader.getTag(tagText, "name");
			graphUnits=XMLReader.getTag(tagText, "units");
	
			graphVarChoice.removeAll();
			for(int j=0;j<graphVarTags.length;j++){
				String s=XMLReader.getTag(graphVarTags[j], "name");
				graphVarChoice.add(s);
			}
		}catch(Exception e){
			Utils.showMessage("Error reading "+varFile+"<br/>"+e);
		}
	}

	public static void readMapVarDataXML(String varFile){
/// Read map variable xml file
		try{
		  String XMLContent=GenieFileReader.getXMLContents(varFile);
			mapVarTags=XMLReader.getTags(XMLContent+"<end></end>", "variable"); //!? end tag needed to fix tag reading error !?
System.out.println("mapVar tags content:");
			for(int i=0;i<mapVarTags.length;i++)System.out.println(mapVarTags[i]);
			mapVarIdx=0;
			setMapFile(mapVarTags[0]);
			readInitialMapVal(0);
	
			mapVarChoice.removeAll();
			for(int j=0;j<mapVarTags.length;j++){
				String s=XMLReader.getTag(mapVarTags[j], "name");
				mapVarChoice.add(s);
			}
		}catch(Exception e){
			Utils.showMessage("Error reading "+varFile+"<br/>"+e);
		}
	}

	public static void setMapFile(String tagText){
		mapTitle=XMLReader.getTag(tagText, "name");
		mapUnits=XMLReader.getTag(tagText, "units");
		String tag=XMLReader.getTag(tagText, "offset");
		offset=0.;
		try{
			offset=Double.valueOf(tag).doubleValue();
		}catch(Exception e){System.out.println(mapTitle+": No offset.");}
		tag=XMLReader.getTag(tagText, "scale");
		scale=1.;
		try{
			scale=Double.valueOf(tag).doubleValue();
		}catch(Exception e){System.out.println(mapTitle+": No rescale.");}
		mapFile=Utils.removeWhitespace(XMLReader.getTag(tagText, "file"));
		initialMapFile=Utils.removeWhitespace(XMLReader.getTag(tagText, "initialFile"));

		nKLayers=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "nKLayers"));
		kLayer=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "kLayer"))-1;
		nLLayers=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "nLLayers"));
		lLayer=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "lLayer"))-1;
		String dirS=XMLReader.getAttribute(tagText, "file", "dirn");
		if(dirS.equals("NSEW"))dirn=0;
		else if(dirS.equals("NSWE"))dirn=1;
		else if(dirS.equals("SNEW"))dirn=2;
		else if(dirS.equals("SNWE"))dirn=3;
		else if(dirS.equals("EWNS"))dirn=4;
		else if(dirS.equals("EWSN"))dirn=5;
		else if(dirS.equals("WENS"))dirn=6;
		else if(dirS.equals("WESN"))dirn=7;
		String typeS=XMLReader.getAttribute(tagText, "variable", "type");
		isOcean3D=false;
		if(XMLReader.getAttribute(tagText, "variable", "type").indexOf("ocean")>=0){
			isOcean=true;
			isAll=isLand=false;
			if(XMLReader.getAttribute(tagText, "variable", "type").indexOf("3D")>=0)isOcean3D=true;
			else depth=7;
		}else if(XMLReader.getAttribute(tagText, "variable", "type").equals("land")){
			isLand=true;
			isOcean=isAll=false;
		}else if(typeS.equals("atmosphere")||typeS.equals("all")){
			isAll=true;
			isOcean=isLand=false;
		}
		if(isOcean3D){
			ocean3DPanel.setVisible(true);
			depth=sDepth;
			kLayer=depth;
		}else ocean3DPanel.setVisible(false);
		nInitialMaps=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "nInitialMaps"));
		nOuterLayers=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "nOuterLayers"));
		if(nOuterLayers>0)outerLayer=Integer.valueOf(XMLReader.getAttribute(tagText, "file", "outerLayer"))-1;
		else outerLayer=0;
	}

	void readInitialData(){
		nInitialMaps=0;
		nOuterLayers=1;
		outerLayer=0;
		nKLayers=8;
		kLayer=0;
		nLLayers=1;
		lLayer=0;
		dirn=0;
		scale=1.;
		offset=0.;
		for(int k=0;k<8;k++){
			kLayer=k;
			if(readMapFile("inputData/seaBed.dat", false).equals("OK")){
				for(int i=0;i<iMax;i++){
					for(int j=0;j<jMax;j++){
						seaBed[i][j][k]=(mapVal[j][i]==1?true:false);
					}
				}
			}
		}
		nKLayers=1;
		kLayer=1;
		nLLayers=1;
		lLayer=1;
		scale=1.;
		offset=0.;
		if(readMapFile("inputData/output.land", false).equals("OK")){
			for(int i=0;i<iMax;i++){
				for(int j=0;j<jMax;j++){
					land[i][j]=(mapVal[j][i]!=0?true:false);
				}
			}
		}
	}
	
	public static String getSpnOrRstFile(String filename){
    File f=new File(filename);
    if(f.exists())return filename;
    if(filename.indexOf("rst.")>=0)filename=filename.replace("rst.","spn.");
    else filename=filename.replace("spn.","rst.");
    File f2=new File(filename);
    if(f2.exists())return filename;
    else return "_NF_";
	}
	
	public static double getCompareTime(String dir){
    compareTime=0;
    String fn=dir+File.separator+"embm/rst.airt";
    fn=getSpnOrRstFile(fn);
    if(fn.equals("_NF_")){
      fn=dir+File.separator+"ents/rst.slandt";
      fn=getSpnOrRstFile(fn);
    }
		String fc=GenieFileReader.getContents(fn);
		int en=fc.lastIndexOf("\n");
		if(en>1)fc=fc.substring(0,en-2);
		int st=fc.lastIndexOf("\n");
		if(st>1)fc=fc.substring(st);
		en=fc.indexOf("  ",4);
    try{
      compareTime=Double.valueOf(fc.substring(0,en)).doubleValue();
    }catch(Exception e){
System.out.println("Error getting last time value from "+fn);
    }
    return compareTime;
	}

	public static void readInitialMapVal(int mapVarI){
    getCompareTime(genieInitialDir);
		readMapFile(genieInitialDir+File.separator+initialMapFile, false);
		for(int i=0;i<iMax;i++)
			for(int j=0;j<jMax;j++)initialMapVal[i][j]=mapVal[i][j];
	}

	public static String readMapFile(String mapFile, boolean change){
		reading=true;
		boolean fileReadOK=true, readFile;
		String fn=getSpnOrRstFile(mapFile);
		if(fn.equals("_NF_")){
System.out.println(mapFile+" & other spn/rst file not found");
      return "Failed to open "+mapFile;
    }
		mapFile=fn;		
		String fileContent=GenieFileReader.getContentsNoBreak(mapFile);
		if(fileContent.equals("")){
			System.out.println("Failed to open "+mapFile);
			reading=false;
			return "Failed to open "+mapFile;
		}else{
  		double num;
			String dataS="";
			StringTokenizer dataST=new StringTokenizer(fileContent, " ");
      mapMin=1e100;
			mapMax=-1e100;
			int i=0, j=0, k=0, l=0, item=0, i0=0, j0=0;
			try{
//if(isOcean3D)kLayer=sDepth;
				for(int ii=0;ii<nInitialMaps;ii++){
		      for (i = iMax1; i >=0; i--){
		        for (j = 0; j <jMax; j++){
							dataS=dataST.nextToken();
						}
					}
				}
        for(int ii=0;ii<nOuterLayers;ii++){
		      for(i=iMax1; i>=0; i--){
		        for(j=0; j<jMax; j++){
		          for(k=0;k<nKLayers;k++){
		            for(l=0;l<nLLayers;l++){
									dataS=dataST.nextToken();
									if(l==lLayer){
										if(k==kLayer){
											if(ii==outerLayer){
												switch(dirn){
													case 0: i0=i; j0=j; break;
													case 1: i0=iMax1-i; j0=j;  break;
													case 2: i0=i; j0=jMax1-j; break;
													case 3: i0=iMax1-i; j0=jMax1-j; break;
													case 4: i0=j; j0=i; break;
													case 5: i0=jMax1-j; j0=i; break;
													case 6: i0=j; j0=iMax1-i; break;
													case 7: i0=jMax1-j; j0=iMax1-i; break;
												}
												if((isLand&&land[j0][i0]) || (isAll) || (isOcean&&!land[j0][i0])){
													num=scale*Double.valueOf(dataS).doubleValue()+offset;
													if(change)num-=initialMapVal[j0][i0];
													mapVal[j0][i0]=num;
					                if(num<mapMin) mapMin=num;
					                if(num>mapMax) mapMax=num;
												}
											}
										}
									}
								}
							}
						}
					}
				}
//System.out.println("nis:"+nInitialMaps+" nos:"+nOuterLayers+" mapVal[0,0]:"+mapVal[0][0]+" mapMin:"+mapMin+" mapMax:"+mapMax);
			}catch(Exception e){
				System.out.print("Possible error reading number from "+mapFile+" "+e);
				reading=false;
				return "Possible error reading number from "+mapFile+" "+e;
			}
      if(mapMin==mapMax){
        mapMin-=1;
        mapMax+=1;
      }
		}
		reading=false;
		return "OK";
	}

////////////
	void drawMap(String s){
		if(mapIm==null || keyIm==null){
			mapIm=createImage(mapW, mapH);
			keyIm=createImage(keyW, mapH);
			return;
		}
		if(s.equals("OK"))drawMap();
		else{
			Graphics mapG=mapIm.getGraphics();
			mapG.setColor(graphPlotter.bgCol);
			mapG.fillRect(0,0, mapW,mapH);
			mapG.setColor(Color.black);
			int l=s.length();
			for(int i=0;i<(int)(l/50)+1;i++){
				int e=(i+1)*50;
				if(e>l)e=l;
				mapG.drawString(s.substring(i*50, e), 10,50+i*20);
			}
			mapPanel.getGraphics().drawImage(mapIm, keyW+2,26, null);
		}		
	}

	void drawMap(){
		if(mapIm==null || keyIm==null){
			mapIm=createImage(mapW, mapH);
			keyIm=createImage(keyW, mapH);
			return;
		}

		Graphics mapG=mapIm.getGraphics();

		mapG.drawString("Waiting for data", 30, 30);
		mapG.drawImage(mapIm, keyW+2,26, null);
		if(reading)return;

		int i, j, i10, j10, layer=7, colI;
		double scl=nMapColrs/(mapMax-mapMin);
		for (i = iMax1; i >=0; i--){
			for (j = 0; j <jMax; j++){
				if(isAll || (isOcean&&!seaBed[i][j][depth]) || (isLand&&land[i][j])){
					colI=(int)((mapVal[i][j]-mapMin)*scl);
					if(colI<0)colI=0;
					else if(colI>nMapColrs1)colI=nMapColrs1;
					mapColrs[i][j]=rainbowCols[colI];
				}else if(land[i][j] && !isLand){
					mapColrs[i][j]=greenCol;
				}else if(!land[i][j] && !isOcean){
					mapColrs[i][j]=blueCol;
				}else if(seaBed[i][j][depth]){
					mapColrs[i][j]=seaBedCol;
				}
			}
		}

		if(globeCB.isSelected())mapG=globePlotter.plotGlobe(mapG, psi, eps, mapColrs, isAll, mapW, mapH);
		else{
			mapG.setColor(graphPlotter.bgCol);
			mapG.fillRect(0,0, mapW,mapH);
			i10=cellW*iMax1;
			for(i=iMax1; i>=0; i--){
				j10=0;
				for(j=0; j<jMax; j++){
					mapG.setColor(mapColrs[i][j]);
					mapG.fillRect(i10,j10, cellW,cellH);
					j10+=cellH;
				}
				i10-=cellW;
			}
			if(isAll)drawCoastline(mapG);
		}
		mapPanel.getGraphics().drawImage(mapIm, keyW+2,26, null);
		drawKey(keyIm.getGraphics(), mapMin, mapMax);

	}

	void drawCoastline(Graphics g){
		Graphics2D g2=(Graphics2D)g;
		g2.setStroke(bs); //set line width
		g2.setColor(greenCol);
	  for(int i=0;i<nCoastlines;i++){
	    g2.drawLine(clStartX[i]*cellW/10, clStartY[i]*cellH/10, clEndX[i]*cellH/10, clEndY[i]*cellH/10);
	  }
	}
////////////

	void drawKey(Graphics g, double min, double max){
		g.setFont(graphPlotter.f);
	  g.setColor(graphPlotter.bgCol);
	  g.fillRect(0,0, keyW,keyH);
	  g.setColor(Color.black);
		double range = max - min;
		if(max==0)max=1e-10;
		if(min==0)min=-1e-10;
		mapValP10=Utils.power(min, max)-1;
		if(mapValP10==1 || mapValP10==-1)mapValP10=0;
		if(mapValP10>1)mapValP10++;
		mapValP=Math.pow(10., mapValP10);
		String s;
//		if(firemaskP->Visible)s="Fire map value";
//		else 
		s=(changeCB.isSelected()?"Change in: ":"")+mapTitle+" / ";
		if(mapValP10!=0)s+="10<sup>"+(mapValP10!=1?""+mapValP10:" ")+"</sup>";
		s+=" "+mapUnits;
		keyTDP.setText(s);
		
		mapDPs=2;
		if(range/mapValP<0.01)mapDPs=4;
		else if(range/mapValP<0.1)mapDPs=3;
		else if(range/mapValP<1)mapDPs=2;
		else if(range/mapValP>10)mapDPs=0;
		else mapDPs=1;
		mapDPP=Math.pow(10, mapDPs);
		for(int i=0;i<=10;i++){
			String ss=""+(int)((min+i*range/10.)/mapValP*mapDPP)/mapDPP;
		  int w=Utils.getWidth(g, ss);
			g.drawString(ss, keyW-17-w, 180-i*17);
		}
		for (int i = 0; i<200; i++){
			g.setColor(rainbowCols[i]);
			g.fillRect(keyW-15, 180-i*18/20, keyW-15, 2);
		}
		keyPanel.getGraphics().drawImage(keyIm, 0, 0, null);
	}

	void readTimeSeriesFile(String filename){
//    Vector<Double> timeValsV=new Vector<Double>();
    
		boolean fileReadOK=true, readFile;
	
		nReadings=-1;
		double num=0;

		String fileContent=GenieFileReader.getContents(filename);
//genieRunDir+File.separator+timeSFile);
		if(fileContent.equals("")){
			System.out.println("Failed to open "+filename);
			graphPlotter.plotBlank("Failed to open "+filename);
//genieRunDir+File.separator+timeSFile);
			fileReadOK=false;
			return;
		}else{
			int stPos=0;
			for(int i=0;i<nHeaderLines;i++)stPos=fileContent.indexOf("\n", stPos+1);
			fileContent=fileContent.substring(stPos).replaceAll("\n", " ");
			StringTokenizer dataST=new StringTokenizer(fileContent, " ");
			int colN=0, nTokens=dataST.countTokens();
//			timeValsV.clear();
			try{
        double oldV=-1e100;
				for(int i=0;i<nTokens;i++){
					colN=i%nCols;
					if(colN==0)nReadings++;
					String dataS=dataST.nextToken();
					if(colN==0){
//            double v=Double.valueOf(dataS).doubleValue();
						timeVal[nReadings]=Double.valueOf(dataS).doubleValue();
//						timeValsV.add(v);
						if(nReadings>0){
							if(timeVal[nReadings]<timeVal[nReadings-1]){
                                Utils.showMessage("Reading time value error: "+
                                        timeVal[nReadings-1]+" "+
                                        timeVal[nReadings]);
							//	nReadings=0;
							//	return;
							}
						}
					}	else if(colN==colNum){
						graphVal[nReadings]=Double.valueOf(dataS).doubleValue();
					}
				}
//				timeVal=timeValsV.toArray();
			}catch(Exception e){
				System.out.print("Error reading number from "+timeSFile+" "+e);
			}
		}
	  if(nReadings>1){
	    time=timeVal[nReadings-1];
			writeGraphData();
	    drawGraph();
	  }else graphPlotter.plotBlank("Waiting for data from "+genieRunDir+File.separator+timeSFile);
	}

	void writeGraphData(){
		String s="<b>Time = "+time+" years</b><br/>";
		if(time>0)s+=graphTitle+" = "
			+Utils.numS(graphVal[nReadings-1], graphPlotter.graphDPs, false)+" "+graphUnits;
		graphTDP.setText(s);
	}

	void writeMouseData(double longitude, double latitude, boolean showVal){
    String s="";
    double LongG, LatP=90-Math.acos(latitude/90)/M_PI*180;

		int i;
    int j=iMax1-(int)(latitude+89.99)/5;

		if(globeCB.isSelected()){
			LongG=longitude;
			i=((int)(longitude+81+180)/10)%36;
		}else{
			i=(int)((longitude+180)/10);
			LongG=longitude-81;
		}
    if(LongG<-180)LongG+=360;
    if(LongG>180)LongG-=360;

    if(showVal){
      double val;
//      if(firemaskP->Visible)s="Fire map value";
//      else 
			s=(changeCB.isSelected()?"Change in: ":"")+mapTitle;
      s=s+" (at "
          +(int)Math.abs(LatP) + "<degrees/> "+(LatP<0?"S":(LatP==0?"":"N"))
          +", "
          +(int)Math.abs(LongG) + "<degrees/> "
					+(LongG<0?"W":(LongG==0?"":"E"))
          +") ";
      boolean landVal=(isOcean ? false : (isLand ? true : land[i][j] ) );
			val=mapVal[i][j];
      if( (land[i][j]&&!landVal)
        ||(!land[i][j]&&landVal)
//        ||(!land[i][j]&&(seaBed[i][j][layer] && fieldPanel->Visible))
			   ){
        s+=": not applicable";
				mapTDP.setText(s);
      }else	mapTDP.setText(s+"=<br/>"+Utils.numS(mapVal[i][j], mapDPs+1, true)+" "+mapUnits);
			mapTDP.setVisible(true);
    }else mapTDP.setVisible(false);
	}

	class MapML extends MouseAdapter{
		public void mouseExited(MouseEvent e){
			if(e.getX()<0 || e.getX()>=mapW || e.getY()<0 || e.getY()>=mapH)mapTDP.setText("");
			e.consume();
		}
		public void mousePressed(MouseEvent e){
		  Y0=e.getY();
		  X0=e.getX();
		}
		public void mouseReleased(MouseEvent e){
		  eps0=eps;
		  psi0=psi;
		}
	}

	class MapMML extends MouseMotionAdapter{
		public void mouseDragged(MouseEvent e){
			if(globeCB.isSelected()&&!reading){
				psi=psi0+(e.getX()-X0)*0.01;
				eps=eps0+(e.getY()-Y0)*-0.01;
				drawMap();
			}
		}
		public void mouseMoved(MouseEvent e){
			int x=e.getX()-keyW-2, y=e.getY()-26;
			if(x>=0 && x<mapW && y>=0 && y<mapH){
				if(x>mapW)x=mapW;
				if(y>mapH)y=mapH;
	      if(globeCB.isSelected()){// && !firemaskP->Visible){
	        double[] pt=globePlotter.calc(x, y, psi, eps);
	        if(pt[0]!=-999){
		        writeMouseData( (int)(pt[0]), -(int)(pt[1]*.5-90)
							, true);
	        } else writeMouseData(1, 1, false);
	      }else{
	        writeMouseData((x*1./mapW-0.5)*360, (y*1./mapH-0.5)*-180, true);
	      }
	    } else writeMouseData(1, 1, false);
			e.consume();
		}
	}
}

