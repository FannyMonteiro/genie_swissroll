import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;

public class Utils{

	public static final double LOG10=Math.log(10);

	public static char nl=(char)(10);
	public static char qt='"';

/** Static finals for quick switching on/off of antialiased drawing:<p>
			Graphics2D g2;
   		g2.getRenderingHints().put( RenderingHints.KEY_RENDERING,
				RenderingHints.VALUE_RENDER_QUALITY);//SPEED);//QUALITY);

			g2.setRenderingHints(gUtils.aaOn);
*/
	public static final float OPAQUE=(float)1.;

	public static final RenderingHints 
				aaOn=new RenderingHints( RenderingHints.KEY_ANTIALIASING,
					RenderingHints.VALUE_ANTIALIAS_ON),
				aaOff=new RenderingHints( RenderingHints.KEY_ANTIALIASING,
					RenderingHints.VALUE_ANTIALIAS_OFF);

/**
Functions: 

	Strings: removeWhiteSpace(s), getWidth(s)
		, drawCenterString(.., drawRightString(..
		, showMessage(s) text in window
		, showAbout(im, s) image & text in window

	Number display: power(min, max), numS(x, nDPs, showPowerOf10)

	Files: 
		showFile(s) text in window
		, readFile(fileName)
		, writeFile(fileName, s)
		, copyFile(inFile, outFile)
		, copyFile(inFile, dir, outFile)

*/


	public static String removeWhitespace(String s){
		s=s.replaceAll(" ", "");
		s=s.replaceAll("\n", "");
		s=s.replaceAll("\t", "");
		return s;
	}

	public static int power(double min, double max){
		int result=0;
		int p1=0, p2=0;
		if(max!=0)p1=(int)(Math.log(Math.abs(max))/LOG10);
		if(min!=0)p2=(int)(Math.log(Math.abs(min))/LOG10);
		result=(p1>p2 ? p1 : p2);
		if(result==1)result=0;
		return result;
	}

	public static String numS(double x, int nDPs, boolean exp){
		if(x==0)return "0";
		int p=(int)(Math.log(Math.abs(x))/LOG10)-1;
		if(p==1 || p==-1)p=0;
		if(p>1)p++;
		x=x/Math.pow(10., p);
		double dp10=Math.pow(10, nDPs);
		return ""+(int)(x*dp10)/dp10+(!exp||(p==0||p==1||p==-1) ? "":"<times/>10<sup>"+p+"</sup>");
	}

	public static int getWidth(Graphics g, String s){
		return g.getFontMetrics().stringWidth(s);
	}

	public static void drawRightString(Graphics g, String s, int x, int y){
		g.drawString(s, x-getWidth(g, s), y);
	}

	public static void drawCenterString(Graphics g, String s, int x, int y){
		g.drawString(s, x-getWidth(g, s)/2, y);
	}

// Show message window
	static Frame frame;
	static String mS;
	static KButton OKB;
	static TextDrawPanel mTDP;

	public static void showMessage(String s){
	  mS=s;
		if(frame==null){
			frame=new Frame();
			frame.setLayout(null);
			frame.setBackground(new Color(255, 255, 240));
			OKB=new KButton("OK", null, true);
    	OKB.setLocation(230,40);
	    OKB.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent e) {
					OKB.down=true;
					OKB.repaint();
					frame.setVisible(false);
					OKB.down=false;
					OKB.repaint();
				}
			});
			frame.add(OKB, null);
			mTDP=new TextDrawPanel();
			mTDP.setBackground(new Color(255, 255, 240));
			mTDP.setBounds(10,70, 480,300);
			mTDP.setExtraTagsEnabled(false);
			mTDP.setFont("SansSerif", Font.PLAIN, 12);
			frame.add(mTDP, null);
		}
    frame.addWindowListener(new WindowAdapter() { 
			public void windowClosing(WindowEvent e) { frame.setVisible(false); } 
		});
    frame.setTitle("Message");
    frame.setSize(500, 200);
    Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
    frame.setLocation((d.width - frame.getSize().width) / 2, (d.height - frame.getSize().height) / 2);

		mTDP.setText(s);
    frame.setVisible(true);
	}

/// About Window
	public static Frame aboutFrame;
	static String aboutMS;
	public static KButton aboutOKB;
	public static TextDrawPanel aboutTDP;
	static Image aboutIm;

	public static void showAbout(String im, String s){
		aboutIm=new ImageIcon(ClassLoader.getSystemResource(im)).getImage();
		int imH=0, imW=500;
		if(aboutIm!=null){
			imH=aboutIm.getHeight(null);
			int w=aboutIm.getWidth(null)+20;
			if(w>imW)imW=w;
		}
	  aboutMS=s;
		if(aboutFrame==null){
			aboutFrame=new Frame(){
				public void paint(Graphics g){
					if(aboutIm!=null)g.drawImage(aboutIm, 10, 30, null);
					super.paint(g);
				}
			};
			aboutFrame.setLayout(null);

			aboutFrame.setBackground(new Color(255, 255, 240));
			aboutTDP=new TextDrawPanel();
			aboutTDP.setBackground(new Color(255, 255, 240));
			aboutTDP.setExtraTagsEnabled(false);
			aboutTDP.setFont("SansSerif", Font.PLAIN, 12);
			aboutFrame.add(aboutTDP, null);

	    aboutFrame.addWindowListener(new WindowAdapter() { 
				public void windowClosing(WindowEvent e) { aboutFrame.setVisible(false); } 
			});
	    aboutFrame.setTitle("About");
	    aboutFrame.setSize(imW, imH+320);
			aboutTDP.setBounds(10,imH+35, imW-20,245);
	    Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
	    aboutFrame.setLocation((d.width - aboutFrame.getSize().width) / 2, (d.height - aboutFrame.getSize().height) / 2);

	    aboutFrame.setVisible(true);

			aboutOKB=new KButton("OK", null, true);
	    aboutOKB.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent e) {
					aboutOKB.down=true;
					aboutFrame.setVisible(false);
					aboutOKB.down=false;
				}
			});
	   	aboutOKB.setLocation(imW/2-20,imH+283);
			aboutFrame.add(aboutOKB, null);

		}

    aboutFrame.setVisible(true);
		aboutTDP.setText(s);
		aboutOKB.revalidate();
  	aboutOKB.repaint();
		aboutOKB.revalidate();
		aboutOKB.requestFocus(true);
	}

// Show file contents in window
	static Frame fFrame;
	static String fS, fFile;
	static TextArea ta;
//	static JButton refreshB;

	public static void showFile(String file){
		fFile=file;
	  fS=readFile(file);
		if(fFrame==null){
			fFrame=new Frame();
			fFrame.setLayout(new BorderLayout());
			fFrame.setBackground(new Color(255, 255, 240));
    	fFrame.setSize(500,500);
			ta=new TextArea(fS);
			ta.setEditable(false);
			fFrame.add(ta, BorderLayout.CENTER);
			JPanel jp=new JPanel();
			jp.setLayout(new BorderLayout());
			JButton refreshB=new JButton("Refresh");
	    refreshB.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent e) {
					fS=readFile(fFile);
					ta.setText(fS);
				}
	    });
			jp.add(refreshB, BorderLayout.EAST);
			JButton hideB=new JButton("Hide");
	    hideB.addActionListener(new ActionListener() {
		    public void actionPerformed(ActionEvent e) {
					fFrame.setVisible(false);
				}
	    });
			jp.add(hideB, BorderLayout.WEST);
			fFrame.add(jp, BorderLayout.NORTH);
		}
    fFrame.addWindowListener(new WindowAdapter() { 
			public void windowClosing(WindowEvent e) { fFrame.setVisible(false); } 
		});
    fFrame.setTitle("File: "+file);
    fFrame.setSize(500,500);
    Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
    fFrame.setLocation((d.width - fFrame.getSize().width) / 2, (d.height - fFrame.getSize().height) / 2);
		ta.setText(fS);
    fFrame.setVisible(true);
	}
	
	public static void centreFrame(Frame f){
    Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
    f.setLocation((d.width - f.getSize().width) / 2, (d.height - f.getSize().height) / 2);
	}
/////////

	public static BufferedReader getFileBR(String inFile){
		BufferedReader dataBR=null;
	  try{
			dataBR=new BufferedReader(new FileReader(inFile));
			System.out.println(inFile+" found");
			return dataBR;
		}catch(Exception e){
			System.out.println(inFile+" not found");
			return null;
		}
	}

  public static String shortenName(String fpS, int len){
    int l=fpS.length(), l2=l;
    while(l>len){
      fpS=fpS.replace(fpS.substring(l/2-8,l/2+8),"...");
      l=fpS.length();
    }
    return fpS;
  }

	public static String readFile(String inFile){
		BufferedReader dataBR=getFileBR(inFile);
		if(dataBR==null) return "";
		char[] data={'c'};
		int dataLength=0;
		String s="", result="";
		StringBuffer sb = new StringBuffer(); 
	  try{
			while(s!=null){
				s=dataBR.readLine();
				if(s!=null)sb.append(s+"\n");
			}
		}catch(Exception e){System.out.println("Error reading "+inFile+": "+dataLength);}
		return sb.toString();
	}

	public static void writeFile(String outFile, String s){
		try{
			FileOutputStream fos=new FileOutputStream(outFile);
			fos.write(s.getBytes());
			fos.close();
		}catch(Exception e){System.out.println("Error wrinting to "+outFile);}
	}

	public static boolean copyFile(String inFile, String outFile){
		File file=new File(inFile);
		if(file.exists()){
			try{
				FileInputStream stream = new FileInputStream(file);
				byte b[] = new byte[stream.available()];
				stream.read(b);
				stream.close();
				File newFile=new File(outFile);
				FileOutputStream output = new FileOutputStream(newFile);
				output.write(b);
				output.flush();
				output.close();
				return true;
			}catch(Exception e){
				System.out.println("Error copying "+inFile+" to "+outFile);
				return false;
			}
		}
		return false;
	}

	public static boolean copyFile(String inFile, String dir, String outFile){
		File file=new File(inFile);
		if(file.exists()){
			try{
				FileInputStream stream = new FileInputStream(file);
				byte b[] = new byte[stream.available()];
				stream.read(b);
				stream.close();
				File f=new File(dir);
				f.mkdir();
				File newFile=new File(outFile);
				FileOutputStream output = new FileOutputStream(newFile);
				output.write(b);
				output.flush();
				output.close();
				return true;
			}catch(Exception e){
				System.out.println("Error copying "+inFile+" to "+outFile);
				return false;
			}
		}
		return false;
	}

}

