import java.awt.*;
import java.awt.event.*;

class MathTest extends Frame{

	String 
		fS="<i>f</i>"
		,rS="<i>r</i>"
		,xS="<i><TRoman>x<SSerif></i>"
		,iS="<i>i</i>"
	;

	TextArea intx=new TextArea("", 0, 0, TextArea.SCROLLBARS_NONE)
	;

	boolean first=true;

	TextDrawPanel tdp=new TextDrawPanel(){
		public void extraTags(){
			TextDraw td=getTextDraw();
			if(td.tag.equals("extratageg/")){
				td.g.drawString("extraTagEg", td.x, td.y);
			}
		}
	};

	public void init(){
		setLayout(null);

		tdp.setBounds(10,140, getWidth()-20, getHeight()-150);
		tdp.setExtraTagsEnabled(true);
		add(tdp, null);

		intx.setBounds(50, 20, 700, 54);
		intx.setText("<br/><br/>abcd ="
			+"<integral>"
			+"<from/><minus/><infinity/> "
			+"<to/>"+xS+"<pow>2x</pow></integral>"
			+fS+"("+rS+") d"+rS+"  "
			+"<fraction>dy<over/>dx</fraction>"
			+"<smaller>"
			+"<sqrt>g<fraction>"
			+"<up><integral>"
			+"<from/>0<to/><infinity/></integral>"
			+"z dz"
			+"<sum><from/>"+iS+"=1 <to/> <infinity/> </sum>"
			+rS+"<sub>"+iS+"</sub></up>"
			+"<over/> <down> <smaller>"+xS+"<pow>3</pow> "
			+"<fraction>  a    <over/>a + b</fraction></smaller></down></fraction>"
			+"abcdefg</sqrt><br/>abc"
		);
		add(intx, null);

		Button rep=new Button("repaint");
		rep.setBounds(0,20,50, 20);
		rep.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(ActionEvent e) {
				tdp.setBounds(10,140, getWidth()-20, getHeight()-150);
				p.setBounds(0,0,getWidth(),getHeight());
				p.repaint();
		}});
		add(rep,null);

		p.setBounds(0,0,getWidth(),getHeight());
		add(p,null);
	}

	Panel p=new Panel(){
		public void paint(Graphics g){
			int sz=(getHeight()-100)/10;
			tdp.setBounds(10,140, getWidth(), getHeight()-140);
			tdp.setFont("SansSerif", Font.PLAIN, sz);
			tdp.setText( intx.getText() );
		}
	};

	public void paint(Graphics g){
		super.paint(g);
		p.repaint();
	}
	
	public static void main(String[] args) {
		MathTest mt=new MathTest();
    mt.addWindowListener(new WindowAdapter() { 
			public void windowClosing(WindowEvent e) { System.exit(0); } 
		});
		mt.setSize(800,500);
		mt.init();
		mt.setVisible(true);
	}
}

