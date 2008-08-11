jar: dirs copy unzip clean pack sign access

dirs:
	rm -rf /tmp/webstart
	mkdir /tmp/webstart
	
copy:	
	cp ../lib/mtj.jar /tmp/webstart/
	cp ../lib/swingintegration0.0.2.jar /tmp/webstart/
	cp ../lib/junit.jar /tmp/webstart/
	cp -rf ../bin/* /tmp/webstart/
	cp -rf ../../CurvatureGeometry/bin/* /tmp/webstart/
	cp -rf ../../Halfedge/bin/* /tmp/webstart/
	cp -rf ../../jreality_halfedge/bin/* /tmp/webstart/
	
unzip:	
	for jar in /tmp/webstart/*.jar;do unzip -o $$jar -d /tmp/webstart; chmod -R +rwx /tmp/webstart/org; done
	
clean:
	rm -f /tmp/webstart/*.jar
	rm -f /tmp/webstart/*.c	
	rm -rf /tmp/webstart/META-INF
	rm -rf /tmp/webstart/CVS
	rm -f /tmp/webstart/LICENSE
	rm -f /tmp/webstart/Makefile
	
pack:
	cd /tmp/webstart; jar cf conformalUnwrapWebstart.jar * 

sign:
	cd /tmp/webstart; jarsigner -storepass hallo123 conformalUnwrapWebstart.jar sechel

access:
	cp -vf /tmp/webstart/conformalUnwrapWebstart.jar ~sechel/WWW/webstart/
	chmod og+r ~sechel/WWW/webstart/conformalUnwrapWebstart.jar
	cp -vf *.jnlp ~sechel/WWW/webstart/
	chmod og+r ~sechel/WWW/webstart/*.jnlp
	rm -rf /tmp/webstart
		
