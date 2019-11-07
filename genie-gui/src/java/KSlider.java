import java.awt.event.*;
import java.awt.*;
import java.net.*;
import java.applet.*;
import java.awt.image.*;
import javax.swing.*;

/**
	Keyboard accessible Slider
*/
public class KSlider extends Button{

	public Image barIm, sliderIm;
	
	int value=0, min=0, x=5, max=100, xMax;

	public void paint(Graphics g){
		super.paint(g);
		g.drawImage(barIm, 0,0,null);
		g.drawImage(sliderIm, x,0,null);
		if(hasFocus()){
			g.setColor(Color.yellow);
			g.drawRect(2, 2, barIm.getWidth(null)-6, barIm.getHeight(null)-6);
			g.drawRect(3, 3, barIm.getWidth(null)-8, barIm.getHeight(null)-8);
		}
	}

	public int getMaximum(){
		return max;
	}

	public int getValue(){
		return value;
	}

	public void setValue(int v){
		value=v;
		x=(v-min)*xMax/(max-min)+5;
		repaint();
	}

	void changeX(){
		if(x<5)x=5;
		else if(x>xMax)x=xMax;
		value=(x-5)*(getMaximum()-min)/xMax+min;
		if(value<min)value=min;
		else if(value>max)value=max;
		repaint();
	}

	public void keyPress(KeyEvent e){
		if(e.getKeyCode()==KeyEvent.VK_TAB){
			e.consume();
			transferFocus();
		}
		if(e.getKeyCode()==KeyEvent.VK_PAGE_UP){
			e.consume();
			x-=(max-min)/10;
			changeX();
		}
		if(e.getKeyCode()==KeyEvent.VK_PAGE_DOWN){
			e.consume();
			x+=(max-min)/10;
			changeX();
		}
		if(e.getKeyCode()==KeyEvent.VK_LEFT){
			e.consume();
			x--;
			changeX();
		}
		if(e.getKeyCode()==KeyEvent.VK_RIGHT){
			e.consume();
			x++;
			changeX();
		}
		if(e.getKeyCode()==KeyEvent.VK_HOME){
			e.consume();
			x=5;
			changeX();
		}
		if(e.getKeyCode()==KeyEvent.VK_END){
			e.consume();
			x=xMax;
			changeX();
		}
	}

	public void mouseDrag(MouseEvent e){
	  x=e.getX();
		changeX();
	}

	public KSlider(int a, int b, Applet ap, boolean isSa, String name){
		super("");
		min=a;
		max=b;
		if(isSa){
			barIm=new ImageIcon("images/"+name+"Bar.png").getImage();
			sliderIm=new ImageIcon("images/"+name+"Slider.png").getImage();
		}else{
			MediaTracker tr=new MediaTracker(ap);
			URL cb=ap.getDocumentBase();
			try {
				barIm = ap.getImage(new URL(cb,"images/"+name+"Bar.png") );
				sliderIm = ap.getImage(new URL(cb,"images/"+name+"Slider.png") );
			}catch (MalformedURLException ex) {}
			tr.addImage(barIm,0);
			tr.addImage(sliderIm,1);
			try {tr.waitForAll();}
			catch (InterruptedException ex){return;}
		}
		setSize(barIm.getWidth(null),barIm.getHeight(null));
		xMax=getWidth()-sliderIm.getWidth(null)-10;
		this.addKeyListener(new KeyListener(){
			public void keyPressed(KeyEvent e){
				keyPress(e);
			}
			public void keyReleased(KeyEvent e){}
			public void keyTyped(KeyEvent e){}
		});
		this.addMouseMotionListener(new java.awt.event.MouseMotionListener() {
			public void mouseMoved(MouseEvent e){	}
			public void mouseDragged(MouseEvent e){
				mouseDrag(e);
			}
			public void mouseReleased(MouseEvent e){	}
		});
	}
}


