import java.awt.*;

public class Colours{

	static final int brightness=200, brightness2=brightness-105;

	public static Color setRainbowCol(double v, double minValue, double maxValue){
		double val=(v-minValue)/(maxValue-minValue);
		double x=1-val;

		int r, g, b;

		if(x<0)x=0;
		else{
		  if(x>1)x=1.;
		}
		if(maxValue<=0)x=(0.5+x/2.);
		else{
		  if(minValue>0)x=x/2.;
		  else{
		    if(maxValue*minValue<=0){
		      double rt=(minValue==0? 1e9 : -maxValue/minValue), k;
		      if(rt<1.){
		        k=0.5-rt*0.5;
		        x=k+x*(1-k);
		      }else{
		        k=0.5+0.5/rt;
		        x=x*k;
		      }
		    }
		  }
		}
		//darken the red
		if (x<.1) {
		    b=0;
		    g=0;
		    r=(int)(brightness2+105*x/.1);
		} else if (x<.4) {
		    b=0;
		    r=brightness;//255;
		    g=(int)((x-.1)/.3*brightness);//255;
		} else if (x<.5) {
		    g=brightness;//255;
		    b=0;
		    r=(int)(brightness-brightness*(x-.4)/.1);
		} else if (x<.9) {
		    r=0;
		    if(x<.6){
		      g=brightness;//255;
		      b=(int)((x-.5)/.1*brightness);//255;
		    }else{
		      g=(int)(brightness-brightness*(x-.6)/.3);
		      b=255;
		    }
		} else {//darken the blue
		    g=0;
		    r=0;
		    b=(int)(255-(x-.9)/.1*150);
		}
		return new Color(r, g, b);
	}
	
}

