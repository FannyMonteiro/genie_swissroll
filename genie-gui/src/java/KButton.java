import java.awt.event.*;
import java.awt.*;
import java.net.*;
import java.applet.*;
import java.awt.image.*;
import javax.swing.*;
import java.lang.ClassLoader;

/**
 Button with Enter key enabled. 
*/

public class KButton extends JButton{
	public KButton() {
		this("");
	}

	public Image upIm, downIm;

	public boolean down=false, enabled=true;
	
	protected void processEvent(AWTEvent e){
		super.processEvent(e);
		repaint();
	}

	public KButton(String text, Applet ap, boolean isSa){
		super(text);
		setSize(text.length()*7+50, 25);
		final Component t=this;
		this.addKeyListener(new KeyAdapter(){
			public void keyPressed(KeyEvent e){
				if (e.getKeyCode() == KeyEvent.VK_ENTER && isEnabled() )
					t.dispatchEvent(new ActionEvent(t, ActionEvent.ACTION_PERFORMED, ""));
			}
		});
		this.addMouseListener(new java.awt.event.MouseListener() {
			public void mousePressed(MouseEvent e){
				down=true;
				repaint();
			}
			public void mouseClicked(MouseEvent e){}
			public void mouseEntered(MouseEvent e){}
			public void mouseExited(MouseEvent e){}
			public void mouseReleased(MouseEvent e){
				down=false;
				repaint();
			}
		});
		Font f=getFont();
		setFont(new Font("SansSerif", Font.BOLD, 11));
		repaint();
 	}

	public KButton(String text){
		super(text);
		final Component t=this;
		this.addKeyListener(new KeyAdapter(){
			public void keyPressed(KeyEvent e){
				if (e.getKeyCode() == KeyEvent.VK_ENTER && isEnabled() )
					t.dispatchEvent(new ActionEvent(t, ActionEvent.ACTION_PERFORMED, ""));
			}
		});
	}
}

