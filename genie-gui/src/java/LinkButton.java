import java.awt.event.*;
import java.awt.*;
import javax.swing.*;

public class LinkButton extends TextDrawPanel{

	String txt;

/*	public void paint(Graphics g){
		super.paint(g);
		g.setColor(this.getBackground());
		g.fillRect(0,0, this.getWidth(),this.getHeight());
		Font f=new Font("SansSerif",Font.BOLD,12);
		g.setFont(f);
		if(hasFocus())g.setColor(Color.blue);
		else 
			g.setColor(Color.black);
		g.drawString(txt, 10,10);
	}//*/

	public LinkButton(String text, boolean isWindows, String osS){
//		super(text);
		txt=text;
		setText(text);
		this.setCursor(new Cursor(Cursor.HAND_CURSOR));
		this.setFocusable(true);
		setSize(text.length()*8+50, 15);
		final Component t=this;
		int p=text.indexOf("http");
		if(p>=0){
			final String pageS=text.substring(p), os=osS;
			final boolean isW=isWindows;
			this.addKeyListener(new KeyAdapter(){
				public void keyReleased(KeyEvent e){
					repaint();
				}
				public void keyPressed(KeyEvent e){
					repaint();
					if (e.getKeyCode() == KeyEvent.VK_ENTER && isEnabled() ){
						try{
							if(isW){
								if (os.equals("Windows 95") || os.equals("Windows 98") || os.equals("Windows ME")) {
									Runtime.getRuntime().exec("start " + pageS);
								} else {
									Runtime.getRuntime().exec("cmd /c start \"name\" \"" + pageS + "\"");
								}
							}	else Runtime.getRuntime().exec("open "+pageS);
						}catch(Exception e1){}
					}
				}
			});
			this.addMouseListener(new java.awt.event.MouseListener() {
				public void mousePressed(MouseEvent e){
					try{
						if(isW){
							if (os.equals("Windows 95") || os.equals("Windows 98") || os.equals("Windows ME")) {
								Runtime.getRuntime().exec("start " + pageS);
							} else {
								Runtime.getRuntime().exec("cmd /c start \"name\" \"" + pageS + "\"");
							}
						}	else Runtime.getRuntime().exec("open "+pageS);
					}catch(Exception e1){}
				}
				public void mouseClicked(MouseEvent e){}
				public void mouseEntered(MouseEvent e){}
				public void mouseExited(MouseEvent e){}
				public void mouseReleased(MouseEvent e){
				}
			});
		}
 	}
}

