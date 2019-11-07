import java.io.File;
import javax.xml.transform.Transformer;
import org.w3c.dom.Document;
import util.xml.*;

public class XSLT
{
	public static void main(String[] args) throws Exception
	{
		if(args.length<3)
		{
			System.err.println("Usage: java XSLT transform.xsl input.xml output.xml");
		}
		Transformer t=XSL.newTransformer(new File(args[0]));
		t.setParameter("user_xml", args[1]);
		t.setParameter("NML_DIR", ".");
		Document d=XSL.transform(t,XML.parse(new File(args[1])));
		XML.save(d, new File(args[2]));
	}
}
