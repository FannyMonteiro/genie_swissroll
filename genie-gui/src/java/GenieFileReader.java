import java.io.*;
import java.util.*;

public class GenieFileReader{

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

	public static String getContents(String inFile){
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

	public static String getXMLContents(String inFile){
		String s=getContents(inFile);
		StringBuffer sb = new StringBuffer(); 
		int offs=0, oldPos=0;
		int pos=s.indexOf("<!--", offs);
		while(pos>=0 && offs>=0){
			sb.append(s.substring(offs, pos));
			offs=s.indexOf("-->", pos);
			pos=s.indexOf("<!--", offs);
		}
		sb.append(s.substring((offs>=0?offs:0)));
		String result=sb.toString();
		result=result.replaceAll("\n", " ").replaceAll(""+Utils.nl, " ").replaceAll(""+(char)(13), " ");
		return result;
	}

	public static String getContentsNoBreak(String inFile){
		BufferedReader dataBR=getFileBR(inFile);
		if(dataBR==null) return "";
		char[] data={'c'};
		int dataLength=0;
		String s="", result="";
		StringBuffer sb = new StringBuffer(); 
	  try{
			while(s!=null){
				s=dataBR.readLine();
				if(s!=null)sb.append(s+" ");
			}
		}catch(Exception e){System.out.println("Error:"+dataLength);}
		return sb.toString();
	}

	public static double[] getDoubles(String inFile){
		Vector<Double> resultV=new Vector<Double>();
		FileInputStream fis=null;
		try{
			fis = new FileInputStream(inFile);
		}catch(Exception e){
			System.out.println(inFile+" found");
			return null;
		}
		try{
			DataInputStream dis = new DataInputStream(fis);
			boolean EOF = false;
			while(!EOF)	{
				try{
					resultV.add(dis.readDouble());
				}catch(Exception e){ EOF = true; }
			}
		}catch(Exception e){ System.out.println("Error reading "+inFile); }
		double result[]=new double[resultV.size()];
		for(int i=0;i<resultV.size();i++){
			result[i]=resultV.get(i);
		}
System.out.println("Nums:"+result[0]+" "+result[1]+" "+result[2]+" "+result[3]+" "+result[4]);
		return result;
	}
}

