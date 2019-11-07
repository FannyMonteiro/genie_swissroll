import java.awt.*;

public class GraphPlotter extends Component{
	/////// Global variables to control graph dimensions
	int gx0=52, gy0=5, gW=285, gH=230, pW, pH;
	double xMin=1, xMax=2
	       , yMin=1e200, yMax=-1e200, gyD, gxD
	       , co2Min=1e200, co2Max=-1e200, gco2D
	       , range
				 , log10=Math.log(10);

	public int graphDPs=2;
	double graphP10=1;

	public int yPower=0;

	public Color bgCol, ltGrayCol;

	public Font f=null;

	Image graphIm;

	Panel gp;

	public String graphTitlePowerS=" ";

	public GraphPlotter(){
	}

	public GraphPlotter(Panel p){
		gp=p;
		pW=p.getWidth();
		pH=p.getHeight();

		graphIm=p.createImage(pW, pH);

		gW=pW-gx0-10;
		gH=pH-gy0-20;
		yMin=1e200;
		yMax=-1e200;
		co2Min=1e200;
		co2Max=-1e200;
		
		log10=Math.log(10);
		
		graphDPs=2;

		graphTitlePowerS=" ";
		
		bgCol=new Color(255, 255, 240);
		ltGrayCol=new Color(180,180,180);
	}

	Graphics graphG;

	//---------------
	int yVal(double y){
	  return (int)(gH+gy0-(y-yMin)*gyD);
	}
	//---------------------------------------------------------------------------
	int co2Val(double y){
	  return (int)(gH+gy0-(y-co2Min)*gco2D);
	}
	//---------------------------------------------------------------------------
	int xVal(double x){
	  return (int)(gx0+(x-xMin)*gxD);
	}
	//---------------------------------------------------------------------------
	int getWidth(Graphics g, String s){
		return g.getFontMetrics().stringWidth(s);
	}

	void drawYlabel(int i, double fy, double p10){
	  double val=i*1./fy;
	  int y=yVal(val);
	  if(y<gy0)return;
	  val=val/p10;
		String valS=""+(int)(val*graphP10+0.5)/graphP10;
	  int w=getWidth(graphG, valS);
	  graphG.setColor(Color.black);
	  graphG.drawString(valS, gx0-3-w, y+4);

	  graphG.drawLine(gx0, y, gx0+4, y);
	  graphG.setColor(ltGrayCol);
	  graphG.drawLine(gx0+4, y, gx0+gW, y);
	}
	//---------------------------------------------------------------------------
	void drawco2label(int i, double fy){
	  double val=i*1./fy;
	  int y=co2Val(val);
	  if(y<gy0)return;
	  int w=getWidth(graphG, Double.toString(val));
	  graphG.drawString(Double.toString(val), gx0+gW+5, y-7);
	  graphG.drawLine(gx0+gW, y, gx0+gW-4, y);
	}
	//---------------------------------------------------------------------------
	public void plotBlank(String s){
		if(graphIm==null){
			graphIm=gp.createImage(pW, pH);
			graphIm.getGraphics().drawString("Waiting for data", 10, 10);
			gp.getGraphics().drawImage(graphIm, 0, 0, null);
			return;
		}
		Graphics g=graphIm.getGraphics();
		g.setColor(bgCol);
		g.fillRect(0,0, pW,pH);
		g.setColor(Color.black);
		int l=s.length();
		for(int i=0;i<(int)(l/50)+1;i++){
			int e=(i+1)*50;
			if(e>l)e=l;
			g.drawString(s.substring(i*50, e), 10,50+i*20);
		}
		gp.getGraphics().drawImage(graphIm, 0, 0, null);
	}

	public void plotGraph(Graphics g, double timeVal[], int xScale, double xMn, double xMx, double yVals[], int nYVals, double yMn, double yMx){
		if(graphIm==null){
			graphIm=gp.createImage(pW, pH);
			return;
		}
		graphG=graphIm.getGraphics();

		if(f==null){
			f = graphG.getFont();
			String fontName = f.getName();
			int fontStyle = f.getStyle();
			f=new Font(fontName, fontStyle, 10);
		}
		graphG.setFont(f);

	  xMin=xMn;
	  xMax=xMx;
	  yMin=yMn;
	  yMax=yMx;
	
	  if(yMin>=1e200 || yMax<=-1e200 || yMax<yMin || xMax<xMin)return;
		if(yMax==yMin){
			yMax+=1;
			yMin-=1;
		}
		if(xMax==xMin){
			xMax+=1;
			xMin-=1;
		}
	
	// Write title
	  range = yMax - yMin;
	  if(yMax==0)yMax=1e-10;
	  if(yMin==0)yMin=-1e-10;

		int p=Utils.power(yMin, yMax)-1;
		if(p==1 || p==-1)p=0;
		if(p>1)p++;

	  double gP10=(Math.pow(10., p));
    graphDPs=2;
    if(range/gP10<0.00001) graphDPs=7;
    else if(range/gP10<0.0001) graphDPs=6;
    else if(range/gP10<0.001) graphDPs=5;
    else if(range/gP10<0.01) graphDPs=4;
    else if(range/gP10<0.1) graphDPs=3;
		graphP10=Math.pow(10., graphDPs);

		if(p!=0)graphTitlePowerS="10<sup>"+(p!=1?""+p:" ")+"</sup>";
		else graphTitlePowerS=" ";

	  gxD=gW/(xMax-xMin);
	  gyD=gH/(yMax-yMin);
	  graphG.setColor(bgCol);
	  graphG.fillRect(0,0,pW, pH);
	  graphG.setColor(Color.black);
	
	// Draw y axis
	  int dy=1, fy=1, p10=1;
	  while(Math.abs(yVal(1*p10)-yVal(0))<gH){
	    if(Math.abs(yVal(1*p10)-yVal(0))<15)dy=2*p10;
	    if(Math.abs(yVal(2*p10)-yVal(0))<15)dy=5*p10;
	    if(Math.abs(yVal(5*p10)-yVal(0))<15)dy=10*p10;
	    p10=p10*10;
	  }
	  p10=1;
	  while(Math.abs(yVal(1./p10)-yVal(0))>gH/5){
	    if(Math.abs(yVal(1./p10)-yVal(0))>40)fy=2*p10;
	    if(Math.abs(yVal(0.5/p10)-yVal(0))>40)fy=5*p10;
	    if(Math.abs(yVal(0.2/p10)-yVal(0))>40)fy=10*p10;
	    p10*=10;
	  }

		yPower=Utils.power(yMin, yMax);
	
	  if(yMin*yMax<=0){
	    for(int i=0;i<=(int)(yMax*fy);i+=dy) drawYlabel(i, fy, gP10);
	    for(int i=dy;i>=(int)(yMin*fy);i-=dy) drawYlabel(i, fy, gP10);
	  }else{
	    if(yMax>0){
	      for(int i=(int)(yMin*fy);i<=(int)(yMax*fy);i+=dy) if(i*1./fy>yMin)drawYlabel(i, fy, gP10);
	    }else{
	      for(int i=(int)(yMax*fy);i>=(int)(yMin*fy);i-=dy) if(i*1./fy<yMax)drawYlabel(i, fy, gP10);
	    }
	  }
	
	  // Draw x axis
	  int w, xw, x, y, dx=1,  fx=1;//x0=(int)xMin,;
	  p10=1;
	  while(xVal(1*p10)-xVal(0)<gW){
	    if(xVal(1*p10)-xVal(0)<40)dx=2*p10;
	    if(xVal(2*p10)-xVal(0)<40)dx=5*p10;
	    if(xVal(5*p10)-xVal(0)<40)dx=10*p10;
	    p10=p10*10;
	  }
	  p10=1;
	  while(Math.abs(xVal(1/p10)-xVal(0))>gW/5){
	    if(Math.abs(xVal(1/p10)-xVal(0))>70)fx=2*p10;
	    if(Math.abs(xVal(0.5/p10)-xVal(0))>70)fx=5*p10;
	    if(Math.abs(xVal(0.2/p10)-xVal(0))>70)fx=10*p10;
	    p10*=10;
	  }
	
	  for(int i=0;i<=(int)(xMax*fx);i+=dx){
	    double val=i*1./fx;
	    if(val>=xMin){
	      x=xVal(val);
				y=yVal(yMin);
				String valS=Double.toString(val);
	  		w=getWidth(graphG, valS);
	      graphG.setColor(Color.black);
	      graphG.drawString(valS, x-w/2, y+10);
	      graphG.drawLine(x, y, x, y-4);
	      graphG.setColor(ltGrayCol);
	      graphG.drawLine(x, y-4, x, y-gH);
	    }
	  }
	  graphG.setColor(Color.black);
	  for(int i=dx;i<=(int)(-xMin*fx);i+=dx){
	    double val=i*1./fx;
	    x=xVal(-val);
			y=yVal(yMin);
			String valS=Double.toString(val);
	  	xw=x-getWidth(graphG, valS)/2;
	    graphG.drawString(valS, xw+5, y+10);
	    graphG.drawLine(xw, y+12, xw+4, y+12);
	    graphG.drawLine(x, y, x, y-4);
	  }
	  graphG.drawLine(gx0, gy0, gx0+gW, gy0);
	  graphG.drawLine(gx0, gy0+gH, gx0+gW, gy0+gH);
	  graphG.drawLine(gx0, gy0, gx0, gy0+gH);
	  graphG.drawLine(gx0+gW, gy0, gx0+gW, gy0+gH);
	
	  //Plot data
	  graphG.setColor(Color.red);;
		int ox=xVal(timeVal[0]*xScale);
		int oy=yVal(yVals[0]);
		Graphics2D g2=(Graphics2D)graphG;
 		g2.getRenderingHints().put( RenderingHints.KEY_RENDERING,
			RenderingHints.VALUE_RENDER_QUALITY);//SPEED);
		g2.setRenderingHints(Utils.aaOn);

	  for(int i=1;i<nYVals;i++){
			x=xVal(timeVal[i]*xScale);
			y=yVal(yVals[i]);
	    graphG.drawLine(ox, oy, x, y);
			ox=x;
			oy=y;
		}

		gp.getGraphics().drawImage(graphIm, 0, 0, null);
	}

}

