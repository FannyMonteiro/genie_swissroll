public class DefnParam extends Object{

	public String name;
	public String value;
	public double min, max, interval;
	public String description;
	public String units;

	public DefnParam(String nm, String v, String mn, String mx, String intvl, String ds, String u){
		name=nm;
		value=v.replaceAll("<sep/>", "/");
		min=Double.valueOf(mn);
		max=Double.valueOf(mx);
		interval=Double.valueOf(intvl);
		while(ds.indexOf("  ")>=0)ds=ds.replaceAll("  ", " ");
		ds=ds.replaceAll("\n", "<br/>");
		description=ds;
		units=u;
 	}

}

