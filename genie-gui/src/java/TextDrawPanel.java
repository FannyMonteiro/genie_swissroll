import java.awt.datatransfer.*;
import java.awt.event.*;
import java.awt.*;

public class TextDrawPanel extends Panel{
	boolean first=true;
	String text="";
	int w, h;
	Color
		background=//new Color(0,0,0,0)
		Color.white
		, borderColor=background;

	Image buffer;
	Graphics bufferG;

	public void paint(Graphics g){
		super.paint(g);
		if(first){
			MediaTracker t=new MediaTracker(this);
			buffer=createImage(w,h);
			t.addImage(buffer,0);
			try{
			t.waitForAll();}catch(Exception e){}
			bufferG=buffer.getGraphics();
			first=false;
		}
//System.out.println(w+" "+h);
		bufferG.setColor(background);
		bufferG.fillRect(0,0,w, h);
		bufferG.setColor(borderColor);
		bufferG.drawRect(0,0,w-1, h-1);
		bufferG.setColor(Color.black);
		bufferG.setFont(tFont);
//*/
/*		g.setColor(background);
		g.fillRect(0,0,w, h);
		g.setColor(borderColor);
		g.drawRect(0,0,w-1, h-1);
		g.setColor(Color.black);
		g.setFont(tFont);
		td.drawPara(g, text, 0, size, this.getWidth());
//*/

		td.drawPara(bufferG, text, 2, size, this.getWidth()-4);
//		ta.setText(getReadableString());
		ta.setLabel(getReadableString());
		g.drawImage(buffer,0,0,null);
	}

	public String getReadableString(){
		return td.getReadableString();
	}

	public void setTabWidth(int t){
		td.setTabWidth(t);
	}

//	public TextArea ta=new TextArea("", 0, 0, TextArea.SCROLLBARS_NONE);
	public Button ta=new Button();

	public TextDrawPanel(){
		super();
		init();
	}//*/

	public TextDrawPanel(String a){
		super();
		text=a;
		init();
	}//*/

	public KButton copyB=new KButton("Copy text", null, true);

	public void init(){
		setLayout(null);
		ta.setBounds(0,0,0,0);//150,20);
//		ta.addKeyListener(new KL());
		addMouseListener(new ML());
		ta.addFocusListener(new FL());
//		ta.setEditable(false);
//		ta.setText(text);
		ta.setLabel(text);
		add(ta, null);

		copyB.addActionListener(new ActionListener(){
		public void actionPerformed(ActionEvent e){
			Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
//			clipboard.setContents(new StringSelection(ta.getText()), null);
			clipboard.setContents(new StringSelection(ta.getLabel()), null);
			copyB.setVisible(false);
		}});
		copyB.setVisible(false);
		add(copyB, null);//*/

		td.setFont(tFont);
	}

	public void setBackground(Color c){
		background=c;
	}

	public void setBounds(int x, int y, int w1, int h1){
		super.setBounds(x, y, w1, h1);
		w=w1;
		h=h1;
	}

////////////////////////////
	int style=Font.PLAIN, size=12;
	Font tFont=new Font("SansSerif", style, size);
	public TextDraw td=new TextDraw();

	public void setFont(String fn, int st, int sz){
		style=st;
		size=sz;
		tFont=new Font(fn, st, sz);
	}

	public void setText(String s){
		text=s;
		repaint();
	}

	public void setWidth(int wd){
		w=wd;
	}

	public void setExtraTagsEnabled(boolean b){
		td.setExtraTagsEnabled(this, b);
	}

	public TextDraw getTextDraw(){
		return td;
	}

	public void extraTags(){
	}

	public TextDrawPanel(String a, int b, int c, int d){
		super();
		text=a;
		init();
	}

	class KL extends KeyAdapter{
		public void keyPressed(KeyEvent e){
			if(e.getKeyCode()==KeyEvent.VK_TAB){
				e.consume();
				transferFocus();
			}else if(e.getKeyCode()>31)e.consume();
		}
		public void keyReleased(KeyEvent e){}
		public void keyTyped(KeyEvent e){}
	}

	class ML extends MouseAdapter{
		public void mouseExited(MouseEvent e){
			if(e.getX()<0 || e.getX()>=w || e.getY()<0 || e.getY()>=h)
			copyB.setVisible(false);
			e.consume();
		}
		public void mousePressed(MouseEvent e){
			if(copyB.isVisible())copyB.setVisible(false);
			else{
				int x=e.getX(), y=e.getY();
				copyB.setLocation((x<w-80 ? x : w-80), (y<h-20 ? y : h-20));
				copyB.setVisible(true);
			}
		}
	}

	class FL extends FocusAdapter{
		public void focusGained(FocusEvent f){
			borderColor=Color.black;
			repaint();
		}
		public void focusLost(FocusEvent f){
			borderColor=background;
			repaint();
		}
	}
}

