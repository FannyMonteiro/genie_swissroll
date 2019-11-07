import java.io.File;

import javax.xml.transform.Transformer;
import org.w3c.dom.Document;
import javax.xml.transform.*;
import javax.xml.transform.dom.DOMResult;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import util.xml.*;

import java.awt.*;
import java.io.*;
import java.util.*;

public class ParamSetting{

	static final String min0="-1", max0="-1", int0="0";

	public ParamSetting(){}

	public Vector<AdjstblParam> adjustablesV;

	public String runTimeRoot="../../genie";
	public String runTimeOutDir=".";
	String[] varNames;
	String[] varValues;

	public void runXSLT(String defnXMLFile, String jobXMLFile, String folder, Object[] adjustedParams){

	  String jobXML=GenieFileReader.getXMLContents(jobXMLFile);
		jobXML=jobXML.replaceAll("\t", "");
		while(jobXML.indexOf("  ")>=0)jobXML=jobXML.replaceAll("  ", " ");

		String[] tmpS=new String[adjustedParams.length];
		for(int i=0;i<adjustedParams.length;i++){
			AdjstblParam ap=(AdjstblParam)adjustedParams[i];
			String name=ap.name;
System.out.println("Adjusting:"+name+" "+ap.value+" NML:"+ap.namelist);
			int bPos=name.indexOf("(");
			int idx=-1;
			if(bPos>0){
				idx=Integer.valueOf(name.substring(name.indexOf("(")+1, name.length()-1));
				name=name.substring(0,bPos);
			}

			int confPos=jobXML.indexOf("/config>");
			int modelPos=jobXML.indexOf("model name="+Utils.qt+ap.modelName+Utils.qt, confPos);
			int stPos=jobXML.indexOf(Utils.qt+name+Utils.qt, modelPos);
			if(ap.value!=null){
				String s0, s1;
				if(idx==-1){
					s0=jobXML.substring(stPos, jobXML.indexOf("</param>",stPos));
					int ep=s0.lastIndexOf(">");
					String vs=s0.substring(ep+1, s0.length());
					s1=s0.replace(vs, ap.value);
					jobXML=jobXML.replace(s0, s1);
				}else{
					s0=jobXML.substring(stPos);
					int iPos=s0.indexOf("index="+Utils.qt+idx+Utils.qt);
					s0=s0.substring(iPos);
					s0=s0.substring(0,s0.indexOf("</param>"));
					int ep=s0.lastIndexOf(">");
					String vs=s0.substring(ep+1, s0.length());
					s1=s0.replace(vs, ap.value);
					jobXML=jobXML.replace(s0, s1);
				}
System.out.println("\nREPLACE:"+s0+"\n"+s1);
			}

		}
		jobXML=jobXML.replaceAll("</param>", "</param>\n");
		jobXML=jobXML.replaceAll("</model>", "</model>\n\n");
		Utils.writeFile("jobTMP.xml", jobXML);

		File f=new File(folder);
		String[] l=f.list();
		for(int i=0;i<l.length;i++){
System.out.println(l[i]);
			File f2=new File(folder+"/"+l[i]);
			if(l[i].indexOf("data_")==0)f2.delete();
//			File f3=new File(folder+"/"+l[i]+".old");
//			if(l[i].indexOf("data_")==0 && l[i].indexOf(".old")<0)f2.renameTo(f3);
		}

// Process xml with xsl files:
		try{
			Transformer t=XSL.newTransformer(new File("nmlBuilding/merge.xsl"));
System.out.println("New Transformer merged.xsl");

			t.setParameter("user_xml", "jobTMP.xml");
			t.setParameter("RUNTIME_ROOT", runTimeRoot);
			t.setParameter("RUNTIME_OUTDIR", runTimeOutDir);

			Document d=XSL.transform(t,XML.parse(new File(defnXMLFile)));
System.out.println("Transformed merged.xsl & "+defnXMLFile);
			XML.save(d, new File("merged.xml"));
System.out.println("Saved merged.xml");
	
			t=XSL.newTransformer(new File("nmlBuilding/toNml.xsl"));
			t.setParameter("NML_DIR", folder);
System.out.println("Transform toNml.xsl & merged.xml");
			d=XSL.transform(t,XML.parse(new File("merged.xml")));
System.out.println("Transform OK: "+d);
			XML.save(d, new File("tmp.xml"));
		}catch (Exception e){
			System.out.println("ERROR processing xslt: "+e);
			Utils.showMessage("ERROR processing xslt: "+e);
		}
	}

	public void setParams(String defnXMLFile, String jobXMLFile, String folder, Object[] adjustedParams){
		getParams(defnXMLFile, jobXMLFile, folder, adjustedParams);
	}

	public Vector getParams(String defnXMLFile, String jobXMLFile, String folder){
		return getParams(defnXMLFile, jobXMLFile, folder, null);
	}

	public Vector getParams(String defnXMLFile, String jobXMLFile, String folder, Object[] paramsToSet){

		adjustablesV=new Vector<AdjstblParam>();
		Vector<DefnParam> resultV=new Vector<DefnParam>();
		adjustablesV.clear();
		resultV.clear();
	  String jobXML=GenieFileReader.getXMLContents(jobXMLFile);
		String jobParamsT=XMLReader.getTag(jobXML, "parameters")+"<t></t>";
		String controlJobT=XMLReader.getTag(jobParamsT, "control");

		Object[] jCParams=getJobParams(controlJobT);

		String defnXML=GenieFileReader.getXMLContents(defnXMLFile);
		String defnParamsT=XMLReader.getTag(defnXML, "parameters");
		String controlDefnT=XMLReader.getTag(defnParamsT, "control");
		String[] defnCFiles=XMLReader.getTags(controlDefnT, "file");

		String varsT=XMLReader.getTag(defnXML, "vars");
		String[] vars=XMLReader.getTags(varsT, "var");
		varNames=new String[vars.length];
		varValues=new String[vars.length];
		for(int i=0;i<vars.length;i++){
			varNames[i]="<varref>"+XMLReader.getAttribute(vars[i], "var", "name")+"</varref>";
			varValues[i]=XMLReader.getTag(vars[i], "var").replaceAll("<sep/>", "/")
				.replaceAll("<varref>RUNTIME_ROOT</varref>", runTimeRoot);
		}

		StringBuffer flagsSB=new StringBuffer();
		String jConfigT=XMLReader.getTag(jobXML, "config")+"<t></t>";
		String jModels[]=XMLReader.getTags(jConfigT, "model");
		String jModelNames[]=new String[jModels.length];
		for(int i=0;i<jModels.length;i++){
			jModelNames[i]=XMLReader.getAttribute(jModels[i], "model", "name");
		}
		String dConfigT=XMLReader.getTag(defnXML, "config")+"<t></t>";
		String dModels[]=XMLReader.getTags(dConfigT, "model");
		String dModelNames[]=new String[dModels.length];
		boolean[] flagState=new boolean[dModels.length];
		for(int i=0;i<dModels.length;i++){
			dModelNames[i]=XMLReader.getAttribute(dModels[i], "model", "name");
			flagState[i]=false;
			for(int j=0;j<jModelNames.length;j++){
				if(jModelNames[j].equals(dModelNames[i]))flagState[i]=true;
			}
		}
		for(int i=0;i<flagState.length;i++){
			String flagName=XMLReader.getAttribute(dModels[i], "model", "flagname");
			if(flagState[i])flagsSB.append(" "+flagName+"=.true.,"+Utils.nl);
			else flagsSB.append(" "+flagName+"=.false.,"+Utils.nl);
		}

		for(int i=0;i<defnCFiles.length;i++){
			String filename=XMLReader.getAttribute(defnCFiles[i], "file", "name");
			String namelist=XMLReader.getAttribute(defnCFiles[i], "namelist", "name");
			Object[] dCParams=getDefnParams(defnCFiles[i], jCParams, namelist, "", paramsToSet);

			StringBuffer sb=new StringBuffer();
			writeParams(sb, namelist, dCParams, (namelist.equals("GENIE_CONTROL_NML")?flagsSB.toString():""));
//			Utils.writeFile(folder+"/"+filename, sb.toString());
		}

		String[] modelsJobT=XMLReader.getTags(jobParamsT, "model");

		String[] modelsDefnT=XMLReader.getTags(defnParamsT, "model");
		for(int i=0;i<modelsDefnT.length;i++){
			String dModelsName=XMLReader.getAttribute(modelsDefnT[i], "model", "name");

			for(int j=0;j<modelsJobT.length;j++){
				if(XMLReader.getAttribute(modelsJobT[j], "model", "name").equals(dModelsName)){
					Object[] jMParams=getJobParams(modelsJobT[j]);

					String filename=XMLReader.getAttribute(modelsDefnT[i], "file", "name");
					StringBuffer sb=new StringBuffer();

					String nmls[]=XMLReader.getTags(modelsDefnT[i], "namelist");
					for(int k=0;k<nmls.length;k++){
						String nmlsN=XMLReader.getAttribute(nmls[k], "namelist", "name");
						Object[] dMParams=getDefnParams(nmls[k], jMParams, nmlsN, dModelsName, paramsToSet);

						writeParams(sb, nmlsN, dMParams, "");
					}

//					Utils.writeFile(folder+"/"+filename, sb.toString());
				}
			}
		}

		return adjustablesV;
	}

	public Object[] getDefnParams(String inS, Object[] newParam, String nmlst, String mdlNm, Object[] aPs){
		Vector<DefnParam> resultV=new Vector<DefnParam>();

// Add params in paramArrays to result
		String[] paramArrays=XMLReader.getTags(inS, "paramArray");
		for(int j=0;j<paramArrays.length;j++){
			String pName=XMLReader.getAttribute(paramArrays[j], "paramArray", "name");
			String params[]=XMLReader.getTags(paramArrays[j], "param");
			for(int k=0;k<params.length;k++){
				String index=XMLReader.getAttribute(params[k], "param", "index");
				String dt=XMLReader.getAttribute(params[k], "value", "datatype");
				String isS=(dt.equalsIgnoreCase("string")?""+Utils.qt:"");
				String value=isS+XMLReader.getTag(params[k], "value")+isS;
				String ds=XMLReader.getTag(params[k], "description");
				String u=XMLReader.getTag(params[k], "units");
				String min=XMLReader.getTag(params[k], "min");
				String max=XMLReader.getTag(params[k], "max");
				String intv=XMLReader.getTag(params[k], "interval");
				if(min.equals("__NOTAG__"))min=min0;
				if(max.equals("__NOTAG__"))max=max0;
				if(intv.equals("__NOTAG__"))intv=int0;

				DefnParam dp=new DefnParam(pName+"("+index+")", value, min, max, intv, ds, u);

				for(int l=0;l<newParam.length;l++){
					JobParam jp=(JobParam)newParam[l];
					if( dp.name.equals( jp.name ) ){
						dp.value=jp.value;
						if( jp.exposed ){
							if(aPs!=null){
								for(int m=0;m<aPs.length;m++){
									AdjstblParam aPSet=(AdjstblParam)aPs[m];
									if(aPSet.name.equals(dp.name) && aPSet.namelist.equals(nmlst)){
										dp.value=aPSet.value;
										if(aPSet.interval!=0){
											dp.min=aPSet.min;
											dp.max=aPSet.max;
											dp.interval=aPSet.interval;
										}
									}
								}
							}else{
								if(jp.interval!=0){
									dp.min=jp.min;
									dp.max=jp.max;
									dp.interval=jp.interval;
								}
								AdjstblParam ap=new AdjstblParam(dp.name, dp.value, ""+dp.min, ""+dp.max, ""+dp.interval, dp.description, dp.units, nmlst, mdlNm);
								adjustablesV.add(ap);
							}
						}
					}
				}

				resultV.add(dp);
			}
		}

// Add params to result
		String paramSrc=XMLReader.ignoreXMLTag(inS, "paramArray");
		String[] params=XMLReader.getTags(paramSrc, "param");
		for(int j=0;j<params.length;j++){
			String pName=XMLReader.getAttribute(params[j], "param", "name");
			String dt=XMLReader.getAttribute(params[j], "value", "datatype");
			String isS=(dt.equalsIgnoreCase("string")?""+Utils.qt:"");
			String value=isS+XMLReader.getTag(params[j], "value")+isS;
			String ds=XMLReader.getTag(params[j], "description");
			String u=XMLReader.getTag(params[j], "units");
			String min=XMLReader.getTag(params[j], "min");
			String max=XMLReader.getTag(params[j], "max");
			String intv=XMLReader.getTag(params[j], "interval");
			if(min.equals("__NOTAG__"))min=min0;
			if(max.equals("__NOTAG__"))max=max0;
			if(intv.equals("__NOTAG__"))intv=int0;

			DefnParam dp=new DefnParam(pName, value, min, max, intv, ds, u);

			for(int k=0;k<newParam.length;k++){
				JobParam jp=(JobParam)newParam[k];
				if( dp.name.equals( jp.name ) ){
					dp.value=jp.value;
					if( jp.exposed ){
						if(aPs!=null){
							for(int m=0;m<aPs.length;m++){
								AdjstblParam aPSet=(AdjstblParam)aPs[m];
								if(aPSet.name.equals(dp.name) && aPSet.namelist.equals(nmlst)){
									dp.value=aPSet.value;
									if(aPSet.interval!=0){
										dp.min=aPSet.min;
										dp.max=aPSet.max;
										dp.interval=aPSet.interval;
									}
								}
							}
						}else{
							if(jp.interval!=0){
								dp.min=jp.min;
								dp.max=jp.max;
								dp.interval=jp.interval;
							}
							AdjstblParam ap=new AdjstblParam(dp.name, dp.value, ""+dp.min, ""+dp.max, ""+dp.interval, dp.description, dp.units, nmlst, mdlNm);
							adjustablesV.add(ap);
						}
					}
				}
			}

			resultV.add(dp);
		}

		return resultV.toArray();
	}

	public Object[] getJobParams(String inS){
		Vector<JobParam> resultV=new Vector<JobParam>();

// Add params in paramArrays to result
		String[] paramArrays=XMLReader.getTags(inS, "paramArray");
		for(int j=0;j<paramArrays.length;j++){
			String pName=XMLReader.getAttribute(paramArrays[j], "paramArray", "name");
			String params[]=XMLReader.getTags(paramArrays[j], "param");
			for(int k=0;k<params.length;k++){
				String index=XMLReader.getAttribute(params[k], "param", "index");
				String exposed=XMLReader.getAttribute(params[k], "param", "exposed");
				String value=XMLReader.getTag(params[k], "param");
				String min=XMLReader.getAttribute(params[k], "param", "min");
				String max=XMLReader.getAttribute(params[k], "param", "max");
				String intv=XMLReader.getAttribute(params[k], "param", "interval");
				if(min.equals("__NOATTRIBUTE__"))min=min0;
				if(max.equals("__NOATTRIBUTE__"))max=max0;
				if(intv.equals("__NOATTRIBUTE__"))intv=int0;

				JobParam jp=new JobParam(pName+"("+index+")", value, min, max, intv, exposed);
				resultV.add(jp);
			}
		}

// Add params to result
		String paramSrc=XMLReader.ignoreXMLTag(inS, "paramArray");
		String[] params=XMLReader.getTags(paramSrc, "param");
		for(int j=0;j<params.length;j++){
			String pName=XMLReader.getAttribute(params[j], "param", "name");
			String exposed=XMLReader.getAttribute(params[j], "param", "exposed");
			String value=XMLReader.getTag(params[j], "param");
			String min=XMLReader.getAttribute(params[j], "param", "min");
			String max=XMLReader.getAttribute(params[j], "param", "max");
			String intv=XMLReader.getAttribute(params[j], "param", "interval");
			if(min.equals("__NOATTRIBUTE__"))min=min0;
			if(max.equals("__NOATTRIBUTE__"))max=max0;
			if(intv.equals("__NOATTRIBUTE__"))intv=int0;

			JobParam jp=new JobParam(pName, value, min, max, intv, exposed);
			resultV.add(jp);
		}

		return resultV.toArray();
	}

	public void writeParams(StringBuffer sb, String namelist, Object[] params, String flags){
		sb.append("&"+namelist+Utils.nl+flags);
		for(int i=0;i<params.length;i++){
			DefnParam dp=(DefnParam)params[i];
			for(int j=0;j<varNames.length;j++){
				dp.value=dp.value.replaceAll(varNames[j], varValues[j]);
			}
			dp.value=dp.value.replaceAll("<varref>RUNTIME_ROOT</varref>", runTimeRoot);
			dp.value=dp.value.replaceAll("<varref>RUNTIME_OUTDIR</varref>", runTimeOutDir);
			sb.append(" "+dp.name+"="+dp.value+","+Utils.nl);
		}
		sb.append("&END"+Utils.nl);
	}

}

