import java.awt.*;

public class GlobePlotter{

	public static int jMax=36, iMax=36
    , iMax1=iMax-1, jMax1=jMax-1, cellSz=360/iMax;

	public static final double PI2=2*Math.PI, PI2360=PI2/360. 
		, PIO2=Math.PI/2;
		
	public static double latScl=5/90.;

	public int[] cellX=new int[4], cellY=new int[4];
	public int jMax2=jMax/2;

	public double x, y, z, r, longt, latt, xz[]=new double[2], yz[]=new double[2]
		, ce, se, cp, sp;
	public int globeR, globeR2, globeR5, xc, yc;
	public int xCoords[][]=new int[iMax][jMax+1];
	public int yCoords[][]=new int[iMax][jMax+1];
	public int zCoords[][]=new int[iMax][jMax+1];

	BasicStroke bs;
	Color greenCol;
	int nCoastlines;
	int[] clStartX=new int[1000], clStartY=new int[1000], clEndX=new int[1000], clEndY=new int[1000];

	public GlobePlotter(BasicStroke bsI, Color gc, int nc, int[] clsx, int[] clsy, int[] clex, int[] cley){
		bs=bsI;
		greenCol=gc;
		nCoastlines=nc;
		clStartX=clsx;
		clStartY=clsy;
		clEndX=clex;
		clEndY=cley;
	}

	double pt[]=new double[3];

	double[] rotatePt(int ix, int iy){
	  double longtd=(ix-270)*PI2360
			, lat=PIO2-Math.acos( (jMax-iy/cellSz)*latScl-1 );
	  double r1=globeR*Math.cos(lat)
	    , x1=r1*Math.sin(longtd), y1=globeR*Math.sin(lat)
	    , z1=r1*Math.cos(longtd), zs;
	  pt[0]= x1*cp+z1*sp;
	  zs=-x1*sp+z1*cp;
	  pt[1]= y1*ce+zs*se;
	  pt[2]=-y1*se+zs*ce;
	  return pt;
	}

	public Graphics plotGlobe(Graphics g, double psi, double eps, Color[][] mapCols, boolean isAll, int mapW, int mapH){
	
    iMax=mapCols.length;
    jMax=mapCols[0].length;
    cellSz=mapW/iMax;
    latScl=cellSz/2/90.;

		ce=Math.cos(eps);
		se=Math.sin(eps);
		cp=Math.cos(psi);
		sp=Math.sin(psi);

		xc=mapW/2;
		yc=mapH/2;
		globeR=xc-10;
		globeR2=globeR*globeR;
		globeR5=globeR/5;
		double rF=1;

		g.setColor(Color.black);
		g.fillRect(0,0, mapW,mapH);

		double xyz[]=new double[3], z;

		int i10=-cellSz;
		for (int i=0; i<iMax; i++){
			i10+=cellSz;
			int j10=0;
			for (int j=0; j<=jMax; j++){
				xyz=rotatePt(i10, j10);
				j10+=cellSz;
				x=(int)(xyz[0]);
				y=(int)(xyz[1]);
				z=(int)(xyz[2]);
				if(z>=0){
					rF=(double)(globeR)/Math.sqrt(x*x+y*y);
					x*=rF;
					y*=rF;
				};
				xCoords[i][j]=(int)(x)+xc;
				yCoords[i][j]=(int)(y)+yc;
				zCoords[i][j]=(int)(z);
			}
		}

		for (int i=0; i<iMax; i++){
			int i1=(i+1)%iMax;
			for (int j=0; j<jMax; j++){
				if(zCoords[i][j]<globeR5){
					int j1=(j+1);
					cellX[0]=xCoords[i][j];
					cellX[1]=xCoords[i1][j];
					cellX[2]=xCoords[i1][j1];
					cellX[3]=xCoords[i][j1];
					cellY[0]=yCoords[i][j];
					cellY[1]=yCoords[i1][j];
					cellY[2]=yCoords[i1][j1];
					cellY[3]=yCoords[i][j1];
					g.setColor(mapCols[i][j]);
					g.fillPolygon(cellX, cellY, 4);
				}
			}
		}

		if(isAll){
			Graphics2D g2=(Graphics2D)g;
			g2.setStroke(bs); //set line width
			g2.setColor(greenCol);
			int xp1, yp1, zp1, xp2, yp2, zp2;
		  for(int i=0;i<nCoastlines;i++){
				xyz=rotatePt(clStartX[i], clStartY[i]);
				xp1=(int)(xyz[0])+xc;
				yp1=(int)(xyz[1])+yc;
				zp1=(int)(xyz[2]);

				xyz=rotatePt(clEndX[i], clEndY[i]);
				xp2=(int)(xyz[0])+xc;
				yp2=(int)(xyz[1])+yc;
				zp2=(int)(xyz[2]);

				if(zp1<=25 && zp2<=25)g2.drawLine(xp1, yp1, xp2, yp2);
		  }
		}

		return g;
	}

	public double resultY[]=new double[2];

	public double[] rotateY(double x, double z, double angle){
		resultY[0]=x*Math.cos(angle)+z*Math.sin(angle);
		resultY[1]=z*Math.cos(angle)-x*Math.sin(angle);
	 	return resultY;
	}

	public double resultX[]=new double[2];

	public double[] rotateX(double y, double z, double angle){
		resultX[0]=y*Math.cos(angle)+z*Math.sin(angle);
		resultX[1]=z*Math.cos(angle)-y*Math.sin(angle);
	 	return resultX;
	}

	public double resultC[]=new double[2];

	public double[] calc(int x2, int y2, double psi, double eps){
		double result[]=new double[2];
		double x0=x2-xc, x02=x0*x0;
		double y0=y2-yc, y02=y0*y0;
		if(x02+y02>=globeR2){
			result[0]=-999;
			return result;
		}
		try{
			double z0=-Math.sqrt(globeR*globeR-x02-y02);
			double ce=Math.cos(eps), se=Math.sin(eps)
				, cp=Math.cos(psi), sp=Math.sin(psi);
		  double yy= y0*ce-z0*se
		  	, z1= y0*se+z0*ce
		  	, xx= x0*cp-z1*sp
		  	, z2= x0*sp+z1*cp
				, lat=Math.asin(yy/globeR)
				, r=globeR*Math.cos(lat)
				, lng=Math.asin(xx/r)
			;
			if(z2<0){
				if(lng>0)lng=Math.PI-lng;
				else lng=-Math.PI-lng;
			}
			result[1]=((jMax-(Math.cos(PIO2-lat)+1)/latScl)*cellSz);
			result[0]=lng/PI2360+9;
		}catch(Exception e){}
		return result;
	}

}

