public class JobParam extends Object{

	public String name;
	public String value;
	public double min, max, interval;
	public boolean exposed;

	public JobParam(String nm, String v, String mn, String mx, String intvl, String expsd){
		name=nm;
		value=v.replaceAll("<sep/>", "/");
		min=Double.valueOf(mn);
		max=Double.valueOf(mx);
		interval=Double.valueOf(intvl);
		exposed=expsd.equals("true");
 	}

}

