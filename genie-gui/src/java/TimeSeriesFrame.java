import java.awt.*;
import java.io.*;
import java.util.*;
import java.awt.event.*;
import javax.swing.*;

public class TimeSeriesFrame extends Frame{

	static final String min0="-1", max0="-1", int0="0";

  File file=null;
	double year[]=new double[2000], val[]=new double[2000];
	int nPts=1;
	double vMax=10, vMin=0;
	double yrMax=10, yrMin=0;

	KButton DefB, CanB, OKB, plotB;

	TextArea TA=new TextArea("", 80, 40);
	TextDrawPanel mTDP=new TextDrawPanel();

	Image graphIm;
	Color graphBgCol;

	int graphW=600, graphH=150, gx0=100, gW=graphW-gx0-80
		, gy0=graphH-40, gH=graphH-45;

	Panel graphPanel;

	boolean windows;

	public void paint(Graphics g){
		super.paint(g);
		if(g!=null){
			if(graphIm!=null)	g.drawImage(graphIm, 10, 350, null);
			else graphIm=createImage(graphW, graphH);
		}
	}

	public void setVisible(boolean vis){
		if(vis)super.setVisible(true);
		else super.setVisible(false);
	}
	
	public String genieRunDir;

	public TimeSeriesFrame(Color bgCol){

		String os = System.getProperty("os.name");
		windows=(os.indexOf("Windows")>=0);

    addWindowListener(new WindowAdapter() { 
			public void windowClosing(WindowEvent e) { setVisible(false); }
		});

    setTitle("Time series file: not selected.");

    MenuBar mb=new MenuBar();
    Menu fileMenu=new Menu("File");
    MenuItem saveMI=new MenuItem("Select file...");
		saveMI.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
        JFileChooser chooser = new JFileChooser();
        chooser.setCurrentDirectory(new File(genieRunDir));
        int returnVal = chooser.showOpenDialog(new Frame());
        if(returnVal == JFileChooser.APPROVE_OPTION) {
          file=chooser.getSelectedFile();
          String fpS=Utils.shortenName(file.getPath(), 45);
          setTitle("Time series file: "+fpS);
          String inpS=Utils.readFile(file.getPath());
          inpS=inpS.replace("-START-OF-DATA-"+Utils.nl, "");
          inpS=inpS.replace("0.0 0.0"+Utils.nl, "");
          inpS=inpS.replace("999999.0 0.0"+Utils.nl
            +"-END-OF-DATA-", "");
          inpS=inpS.replace(""+Utils.nl+Utils.nl, "");
          TA.setText(inpS);
          
          calc();
          drawGraph();
        }
			}
		});
		fileMenu.add(saveMI);
		mb.add(fileMenu);
		this.setMenuBar(mb);//*/

    setSize(700, 540);
    Utils.centreFrame(this);
    setVisible(false);

		setLayout(null);
		setResizable(false);

		setBackground(bgCol);
		graphBgCol=bgCol;

		mTDP.setBackground(new Color(255, 255, 240));
		mTDP.setBounds(150,55, 400,42);
		mTDP.setExtraTagsEnabled(false);
		mTDP.setFont("SansSerif", Font.PLAIN, 12);
		mTDP.setText("Enter pairs of numbers separated by spaces "
			+"('year[space]value')."
			+"Add each pair on a new line in the box below. "
//			+"(Use no more than 6 significant figures)"
		);
		add(mTDP, null);

		TA.setBounds(150,100, 400,200);
		add(TA, null);

		plotB=new KButton("Plot", null, true); // images folder must have an image 'plotUp.png'
    plotB.setLocation(350-plotB.getWidth()/2,310);
    plotB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				plotB.down=true;
				plotB.repaint();

				calc();
				drawGraph();

				plotB.down=false;
				plotB.repaint();
			}
    });
		add(plotB, null);

/*		DefB=new KButton("Set default values", null, true); 
    DefB.setLocation(100,500);
    DefB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				DefB.down=true;
				DefB.paint(DefB.getGraphics());

				drawGraph();

				DefB.down=false;
				DefB.paint(DefB.getGraphics());
			}
    });
		add(DefB, null);//*/

		CanB=new KButton("Cancel", null, true); 
    CanB.setLocation(250,500);
    CanB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				CanB.down=true;
				CanB.paint(CanB.getGraphics());
				setVisible(false);
				CanB.down=false;
				CanB.paint(CanB.getGraphics());
			}
    });
		add(CanB, null);

		OKB=new KButton("OK", null, true); 
    OKB.setLocation(420,500);
    OKB.addActionListener(new ActionListener() {
	    public void actionPerformed(ActionEvent e) {
				OKB.down=true;
				OKB.paint(OKB.getGraphics());

				if(calc()<=0){
          return;
        }
				String outS="-START-OF-DATA-"+Utils.nl+(year[1]>0?"0.0      0.0":"")+Utils.nl;
			  for(int i=1;i<nPts+1;i++)outS+=year[i]+"    "+val[i]+Utils.nl;
			  outS+="999999.0 0.0"+Utils.nl
			    +"-END-OF-DATA-"+Utils.nl+Utils.nl;
			  outS=outS.replaceAll(""+Utils.nl+Utils.nl, ""+Utils.nl);
			    
        Utils.writeFile(file.getPath(), outS);
        Utils.writeFile(genieRunDir+"/"+file.getName(), outS);
        setVisible(false);

				OKB.down=false;
				OKB.paint(OKB.getGraphics());
			}
    });
		add(OKB, null);

	}

	Font f;

	public void drawGraph(){
		Graphics g=graphIm.getGraphics();
		if(f==null){
			f = g.getFont();
			String fontName = f.getName();
			int fontStyle = f.getStyle();
			f=new Font(fontName, fontStyle, 10);
		}
		g.setFont(f);
		if(g==null)return;
	  if(yrMax==yrMin)yrMax++;
	  if(vMax==vMin)vMax++;
		g.setColor(graphBgCol);
		g.fillRect(0,0, graphW,graphH);
	  g.setColor(Color.black);
		Utils.drawRightString(g, ""+vMax, gx0-5, gy0-gH+4);
		Utils.drawRightString(g, ""+vMin, gx0-5, gy0+3);
		Utils.drawCenterString(g, ""+yrMin, gx0,gy0+17);
		Utils.drawCenterString(g, ""+yrMax, gx0+gW,gy0+17);
	  g.drawLine(gx0,gy0-gH, gx0-3,gy0-gH);
	  g.drawLine(gx0, gy0, gx0-3,gy0);
	  g.drawLine(gx0, gy0, gx0,gy0+3);
	  g.drawLine(gx0+gW, gy0, gx0+gW,gy0+3);
	  g.drawLine(gx0, gy0-gH, gx0, gy0);
	  g.drawLine(gx0,gy0, graphW, gy0);
		Graphics2D g2=(Graphics2D)g;
 		g2.getRenderingHints().put( RenderingHints.KEY_RENDERING,
			RenderingHints.VALUE_RENDER_QUALITY);//SPEED);
		g2.setRenderingHints(Utils.aaOn);
		g2.setColor(Color.red);
		int ox=gx0
			, oy=(int)(gy0-(0-vMin)/(vMax-vMin)*gH);
	  for(int i=1;i<nPts+1;i++){
			int x=(int)((year[i]-yrMin)/(yrMax-yrMin)*gW+gx0)
				, y=(int)(gy0-(val[i]-vMin)/(vMax-vMin)*gH);
	    g2.drawLine(ox,oy, x,y);
			ox=x;
			oy=y;
	  }
	  g2.drawLine(ox,oy, (int)((yrMax+1100-yrMin)/(yrMax-yrMin)*gW+gx0), gy0);
		repaint();
	}

	public int calc(){
	  double num1, num2;
	  year[0]=0;
	  val[0]=0;
		String txt=TA.getText(), lines[]=new String[10000];
	  nPts=0;
		String nTxt="";
		try{
			for(int i=0;i<txt.length();i++){
				char c=txt.charAt(i);
				if(c<32){
				  if(c==Utils.nl)nTxt+="L";
				}else{
					nTxt+=""+c;
				}
			}
			nTxt+="L";
			while(nTxt.indexOf("LL")>=0)nTxt=nTxt.replaceAll("LL","L");
			int pos=0, nPos=nTxt.indexOf("L");
			while(nPos>0){
				lines[nPts]=nTxt.substring(pos, nPos).replaceAll("\t"," ");
				pos=nPos+1;
				nPos=nTxt.indexOf("L", pos);
				nPts++;
			}
			for(int j=0;j<nPts;j++)System.out.println(lines[j]);
		  vMax=0;
		  vMin=0;
		  yrMax=0;
		  yrMin=0;
		  for(int i=1;i<=nPts;i++){
		    String s=lines[i-1];
				while(s.indexOf("  ")>=0)s=s.replaceAll("  "," ");
				s=s.replaceAll(" ", ",");
				if(s.charAt(0)==',')s=s.substring(1,s.length());
				int cp=s.indexOf(",");
				String s1=s.substring(0, cp);
				num1=Double.valueOf(s1);
				int cp2=s.indexOf(",", cp+1);
				String s2=s.substring(cp+1, (cp2>0?cp2:s.length()));
				num2=Double.valueOf(s2);
		    if((i>1 && num1<=year[i-1]) || (i==1 && num1<0)){
		      Utils.showMessage("Please give a year greater than "+year[i-1]+"\n for pair "+i);
		      return 0;
		    }
		    if(num1>=999999.){
		      Utils.showMessage("Please give a year less than 999999\n for pair "+(i));
		      return 0;
		    }
		/*    if(num2<0){
		      ShowMessage("Please give a value greater than 0\n for line "+(String)(i));
		      return 0;
		    }//*/
		    if(num1<yrMin)yrMin=num1;
		    if(num1>yrMax)yrMax=num1;
		    if(num2<vMin)vMin=num2;
		    if(num2>vMax)vMax=num2;
		    year[i]=num1;
		    val[i]=num2;//*/
		  }
		}catch(Exception e){
			Utils.showMessage("Please type in valid pairs of numbers.<br/>"+e);
			return 0;
		}
	  drawGraph();
	  return nPts;
	}
}

