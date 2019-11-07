import java.awt.*;
import javax.swing.*;

public class TextDraw extends Component{

/**
*/

/**
*/
	public Image textIm;

/**
*/
	public int x, xOffset=0;

/**
*/
	public double lineSpacing=1.;

/**
*/
	public int tabWidth=30;

/**
Draws a paragraph of a text string, s, at (x, y) within a width, w, to
Graphics, g,
returning the y coord at which the text ends.
Formatting codes can be included in s.<p>
E.g. to draw the integral from -infinity to (x^2) of
f(r) dr followed by an image in im<p>
type:<br>
*/

	public String readableString="";

	public String getReadableString(){
		return readableString;
	}

	public int getX(){return x;}

	public int getWidth(Graphics g, String s){
		return g.getFontMetrics().stringWidth(s);
	}

	public void setTabWidth(int w){
		tabWidth=w;
	}

	public Font f;
	public String name, tag, params;
	public int style, size, oldSize, xEnd, h;
	public Graphics g;
	public int fromStart[]=new int[100], nFroms=0
			, byStart[]=new int[100], byEnd[]=new int[100], nBys=0
			, sqXmin[]=new int[100], sqYmin[]=new int[100], sqYmax[]=new int[100]
			, nSqs=0
			,x0, y, w
	;

	boolean doExtraTags=false, doExtraJTags=false;

	public int drawPara(Graphics g1, String s, int x01, int y1, int w1){
		g=g1;
		x=x01;
		x0=x01;
		y=y1;
		w=w1;

		s = " "+s + "<endPara> ";

		f = g.getFont();
		name = f.getName();
		style = f.getStyle();
		size = f.getSize();
		oldSize=size;

		readableString="";
		x=x0+xOffset;
		xOffset=0;
		xEnd = x0 + w;

		h = (int)(g.getFontMetrics().getHeight()*lineSpacing);
		
		boolean draw = true;
		String word, text;
		int i = 0;

		try{
			while (i >= 0){

				int tagStart=s.indexOf('<', i + 1)
					, tagEnd=s.indexOf('>', i + 1)
					, tagSpaceStart=s.indexOf(' ', tagStart+1)
				;

				text=s.substring(i+1, tagStart)+" ";

				int j=0, wordEnd=text.indexOf(' ', j+1);
				while(wordEnd>0){
					wordEnd=text.indexOf(' ', j+1);
					if(wordEnd>0){
						word=text.substring(j, wordEnd);
						int length = getWidth(g, word);
						if (x + length > xEnd){
							y += h;
							x = x0-getWidth(g, " ");
							String wd2=word;
							int l2=getWidth(g, wd2);
							if(l2>w){
                while(l2>w){
                  wd2=wd2.substring(0,wd2.length()-2);
                  l2=getWidth(g, wd2);
                }
                y-=h;
                g.drawString(wd2, x,y);
                y+=h;
                word=word.substring(wd2.length());
							}
						}
						if((style&Font.ITALIC)!=0) g.drawString(word+" ", x, y);
						else g.drawString(word, x, y);
						x += length;
						j=wordEnd;
					}
				}
				readableString+=text;

				if(tagSpaceStart<tagEnd){
					tag=s.substring(tagStart+1, tagSpaceStart).toLowerCase();
					params=s.substring(tagSpaceStart+1, tagEnd).toLowerCase();
				}else{
					tag=s.substring(tagStart+1, tagEnd).toLowerCase();
					params="";
				}

				doTags();
//				if(doExtraJTags)tdjp.extraTags();
//				else 
				if(doExtraTags)tdp.extraTags();

				g.setFont(new Font(name, style, size));
				i = tagEnd;
				if(tag.equals("endpara"))i=-1;
			}
		}catch(Exception e){
			System.out.println(e);
			readableString+="Cannot format remaining text";
			g.drawString("Cannot format remaining text", x, y);
		}//*/

		g.setFont(f);
		return y;
	}

//	TextDrawJPanel tdjp;
	TextDrawPanel tdp;

/*	public void setExtraTagsEnabled(TextDrawJPanel tjp, boolean b){
		tdjp=tjp;
		doExtraJTags=b;
	}//*/

	public void setExtraTagsEnabled(TextDrawPanel tp, boolean b){
		tdp=tp;
		doExtraTags=b;
	}
	
	public void doTags(){
		if(tag.equals("<")){
			g.drawString("<", x, y);
			x+=size>>2;
		}
		if(tag.equals("bar/")){
			readableString+=" bar ";
			g.drawString(""+(char)175, x-(size>>1), y);
			x+=2;
		}
		if(tag.equals("i"))	style |= Font.ITALIC; 
		if(tag.equals("/i"))style &= ~Font.ITALIC;
	
		if(tag.equals("sup")||tag.equals("pow")){
			readableString+=(tag.charAt(0)=='p' ? " to the power " : " superscript ");
			size -= 3; y -= size>>1; 
		}
		if(tag.equals("/sup")||tag.equals("/pow")){
			readableString+=(tag.charAt(1)=='p' ? " end power ":" unsuperscript ");
			y += size>>1; 
			size += 3; 
		}
		if(tag.equals("sub")){
			readableString+=" subscript ";
			size -= 3;
			y += size>>1; 
		}
		if(tag.equals("/sub")){
			readableString+=" end subscript ";
			y -= size>>1; 
			size += 3; 
		}
		if(tag.equals("troman")) name = "TimesRoman";
		if(tag.equals("sserif")) name = "SansSerif";
		if(tag.equals("br/")){
//			readableString+=" \n ";
			readableString+=" new line\n ";
			x = x0;
			y += h;
		}
		if(tag.equals("back/")) x -= (int)(size*.5);
		if(tag.equals("b"))	style |= Font.BOLD;
		if(tag.equals("/b")) style &= ~Font.BOLD;
		if(tag.equals("integral")){
			readableString+=" the integral ";
			x+=size>>2;
			g.setFont( new Font( name, Font.PLAIN, (int)(size*1.3) ) );
			g.drawString(""+(char)8747, x, y+(size>>3));
			g.setFont(new Font(name, style, size));
			x+=size>>1;
		}
		if(tag.equals("sum")){
			readableString+=" the sum ";
			x+=size>>2;
			g.setFont( new Font( name, Font.PLAIN, size ) );
			g.drawString(""+(char)8721, x, y-size/7);
			g.setFont(new Font(name, style, size));
		}

		if(tag.equals("infinity/")){
			readableString+=" infinity ";
			g.setFont(new Font(name, Font.PLAIN, size));
			g.drawString(""+(char)8734, x, y);
			g.setFont(new Font(name, style, size));
		}
		if(tag.equals("minus/")){
			readableString+=" minus ";
			g.drawString(""+(char)8722, x, y);
			x+=(size>>1)+3;
		}
		if(tag.equals("degrees/")){
			readableString+=" degrees ";
			g.drawString("\u00b0", x, y);
			x+=(size>>1)+1;
		}
		if(tag.equals("times/")){
			readableString+=" times ";
			g.drawString("\u00d7", x, y);
			x+=(size>>1)+1;
		}

		if(tag.equals("from/")){
			readableString+=" from ";
			oldSize=size; size=(int)(size*.5); y+=(int)(size*.8); 
			//x-=size; 
			fromStart[nFroms++]=x; 
		}
		if(tag.equals("to/")){
			readableString+=" to ";
			y-=(int)(size*2.6); x=fromStart[nFroms-1];//-size;
			nFroms--;
		}
		if(tag.equals("/integral")||tag.equals("/sum")){
			readableString+=" of ";
			y+=(int)(size*1.8); x+=(int)(1.1*size);
			size=oldSize;
		}
		if(tag.equals("larger")||tag.equals("/smaller")){
			size += 2;
		}
		if(tag.equals("smaller")||tag.equals("/larger")){
			size -= 2;
		}
		if(tag.equals("up")||tag.equals("/down")){
			y-=h/6;
		}
		if(tag.equals("down")||tag.equals("/up")){
			y+=h/6;
		}
		if(tag.equals("draw image")){
			if(textIm!=null){
				g.drawImage(textIm, x, y-h*3/4, null);
				x+=textIm.getWidth(null); 
			}
		}
		if(tag.equals("tab/")){
			readableString+="	";
			x=(x+tabWidth)/tabWidth*tabWidth;
		}
		if(tag.equals("fraction")){
			nBys++;
			readableString+=" fraction ";
			size -= 3; y -= (int)(size*.7);
			byStart[nBys]=x;
			x+=(size>>2);
		}
		if(tag.equals("over/")){
			readableString+=" over ";
			y+=(int)(size*1.4); 
			byEnd[nBys]=x+(size>>2);
			x=byStart[nBys]+(size>>2); 
		}
		if(tag.equals("/fraction")){
			if(x>byEnd[nBys])byEnd[nBys]=x;
			int yy=y-size;
			g.drawLine(byEnd[nBys], yy, byStart[nBys], yy); 
			readableString+=" end fraction ";
			x=byEnd[nBys]+(size>>2);
			y-=(int)(size*.7); size += 3;
			nBys--;
		}
		if(tag.equals("</")){
			readableString+="<";
			g.drawString("<", x, y);
			x+=getWidth(g, "<");
		}
		if(tag.equals("sqrt")){
			readableString+=" square root ";
			nSqs++;
			sqXmin[nSqs]=x+5;
			x+=10;
			sqYmin[nSqs]=y;
			sqYmax[nSqs]=y;
		}
		if(tag.equals("/sqrt")){
			readableString+=" end square root ";
			x+=3;
			int sw=x-sqXmin[nSqs]
				, sh=sqYmax[nSqs]-sqYmin[nSqs]
				, xps[]={sqXmin[nSqs]-10, sqXmin[nSqs]-6, sqXmin[nSqs]-3, sqXmin[nSqs], x}
				, yps[]={sqYmax[nSqs]-10, sqYmax[nSqs]-12, sqYmax[nSqs], sqYmin[nSqs], sqYmin[nSqs]}
			;
			g.drawPolyline(xps, yps, 5);
			nSqs--;
		}
		if(tag.equals("table")){
			nTables++;
			tableX[nTables]=x;
			tableY[nTables]=y-size;
			nRows[nTables]=0;
		}
		if(tag.equals("tr")){
			nRows[nTables]++;
			x=tableX[nTables];
			y+=h;
			rowY[nTables][nRows[nTables]]=y;
		}//*/
		if(tag.equals("/td")){
			nCols[nTables]++;
			colX[nTables][nCols[nTables]]=x;
		}
		if(tag.equals("td")){
			int wp=params.indexOf("width=");
			if(wp>0){
//System.out.println(wp);
			}
			nCols[nTables]++;
			colX[nTables][nCols[nTables]]=x;
		}
		if(nSqs>0){
			if(x<sqXmin[nSqs])sqXmin[nSqs]=x;
			if(y+(size>>2)>sqYmax[nSqs])sqYmax[nSqs]=y+(size>>2);
			if(y-(size)<sqYmin[nSqs])sqYmin[nSqs]=y-(size);
			for(int i=0;i<nSqs;i++){
				sqYmax[i]=sqYmax[nSqs]+(i<<1);
				sqYmin[i]=sqYmin[nSqs]-(i<<1);
			}
		}
	}

	int nTables=0, nRows[]=new int[50], tableX[]=new int[50]
		, tableY[]=new int[50], nCols[]=new int[50]
		, rowY[][]=new int[50][50], colX[][]=new int[50][50]
	;
}

