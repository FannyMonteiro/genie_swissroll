import java.util.*;

public class XMLReader{
	/* 
	getTag(inS, tag) returns inner content of tag within string string inS.
	getTags(inS, tag) returns array of tag content within string inS (includes tag itself).
	getAttribute(inS, tag, att) returns 1st occuring attribute value (inside "'s) of att
		from tag within inS.
	st & fn are start & finish positions within a string.
	Maximum no. tags collected = 120000, but can be changed
	*/

	public static int tagIndex(String inS, String tg, int offs){
		int st=inS.indexOf("<"+tg+" ", offs);
		int st2=inS.indexOf("<"+tg+">", offs);
		if(st>=0){
			if(st2>=0){
				if(st2<st)return st2;
				else return st;
			}else return st;
		}
		if(st2>=0)return st2;
		st2=inS.indexOf("<"+tg+"/", offs);
		if(st2>=0)return st2;
		return -1;
	}

	public static String getTag(String inS, String tag){
		if(inS==null || tag==null)return "__NOTAG__";
		String content=getTagO(inS, tag, 0);
		int st=content.indexOf(">")+1;
		int fn=content.lastIndexOf("<");
		if(st>0 && fn>=0)return content.substring(st, fn);
		else return "__NOTAG__";
	}
	
	public static String getTagO(String inS, String tag, int offst){
		int st=tagIndex(inS, tag, offst);
		int fn1=inS.indexOf("/>", st);
		int fn2=inS.indexOf("<", st+1);
		int f=tagIndex(inS, "/"+tag, offst);
		if((fn1<f || (f<0&&fn1>=0)) && fn1>0 && st>=0 && fn1<fn2){
			return inS.substring(st, fn1+2);
		}
		if(st<0 || f<1)return "__NOTAG__";
		f=(inS.indexOf(">", f))+1;
		if(f<1)return "__NOTAG__";
		String content=inS.substring(st, f);
		int tp=tagIndex(content, tag, offst+1);
int cc=0;
		while(tp>=0 && cc<120000){
			f=tagIndex(inS, "/"+tag, f+1);
			f=(inS.indexOf(">", f))+1;
			if(f<1){
				tp=-1;
				content="__NOTAG__";
			}else{
				content=inS.substring(st, f);
				tp=tagIndex(content, tag, tp+1);
			}
cc++;
		}
		return content;
	}

	public static String[] getTags(String inS, String tag){
		Vector<String> resultV=new Vector<String>();
		if(inS==null || tag==null){
			String[] result0=resultV.toArray(new String[resultV.size()]);
			return result0;
		}
		int i=0;
		int l=inS.length();
		String content=getTagO(inS, tag, 0);
		int offset=tagIndex(inS, tag, 0);//+s.length;

int c=0;
		while(offset>=0 && c<120000){
c++;
			resultV.add(content);
			i++;
			if(offset<0)content="__NOTAG1__";
			else{
				inS=inS.substring(offset+content.length());
				content=getTagO(inS, tag, 0);
				offset=tagIndex(inS, tag, 0);//+1;
			}
		}
		String[] result=resultV.toArray(new String[resultV.size()]);
		return result;
	}

	public static String getAttribute(String inS0, String tg, String attr){
		attr=" "+attr+"=";
		String inS=getTagO(inS0+"<>", tg, 0);
		int i=0;
		String resultS="";
		int fn=inS.indexOf(">");
		if(fn<0)return resultS;
		String tagS=inS.substring(0, fn+1);
		int st=tagS.indexOf(attr);
		if(st<0)return "__NOATTRIBUTE__";
		st=tagS.indexOf("\"", st);
		fn=tagS.indexOf("\"", st+1);
		if(fn<0)fn=tagS.indexOf(">", st);
		if(fn<0)fn=tagS.indexOf("/", st);
		if(st<0 || fn<0)return "__NOATTRIBUTE__";
		resultS=tagS.substring(st+1, fn);
		return resultS;
	}

	public static String ignoreXMLTag(String s, String tag){
		StringBuffer sb = new StringBuffer(); 
		int offs=0, oldPos=0;
		int pos=s.indexOf("<"+tag, offs);
		while(pos>=0 && offs>=0){
			sb.append(s.substring(offs, pos));
			offs=s.indexOf("</"+tag, pos);
			pos=s.indexOf("<"+tag, offs);
		}
		sb.append(s.substring((offs>=0?offs:0)));
		return sb.toString();
	}

}

