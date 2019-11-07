import java.util.*;
import java.io.*;

public class HTMLFiles{

	public HTMLFiles(){}

	String mS[]={"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

	public String runFolder="";

	public void createHTMLFiles(String paramS, String folder){
		String ts=(new Date()).toString();
		String dayS=ts.substring(8,10);
//		ts=ts.substring(3);
		String monthS="";
		for(int i=1;i<=mS.length;i++)if(ts.indexOf(mS[i-1])>=0)monthS=(i<10?"0":"")+i;

		ts="runOn_"+dayS+"_"+monthS+"_"+ts.substring(24, 28)+"_at_"+ts.substring(11, 23);
		runFolder=folder+"/"+ts.replaceAll(":",".").replaceAll(" ","_");
System.out.println("\n\nSaving html to "+runFolder);
		File f=new File(runFolder);
		f.mkdir();
		Utils.writeFile(runFolder+"/parameterValues.html", paramS);
	}

}
