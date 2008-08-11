jar: dirs copy unzip clean pack sign access

dirs:
	rm -rf /tmp/webstart
	mkdir /tmp/webstart
	
copy:	
#	cp ../lib/jReality.jar /tmp/webstart/
	cp ../lib/jSyntax.jar /tmp/webstart/
#	cp ../lib/synthetica.jar /tmp/webstart/
#	cp ../lib/syntheticaBlueSteel.jar /tmp/webstart/
#	cp ../lib/syntheticaGreenDream.jar /tmp/webstart/
#	cp ../lib/syntheticaSilverMoon.jar /tmp/webstart/
	cp ../lib/mtj.jar /tmp/webstart/
#	cp ../lib/jtem-beans.jar /tmp/webstart/
	cp ../lib/numericalMethods.jar /tmp/webstart/
	cp ../lib/vecmath.jar /tmp/webstart/
	cp ../lib/mfc.jar /tmp/webstart/
#	cp ../lib/bsh.jar /tmp/webstart/	
#	cp ../lib/itext-1.4.6.jar /tmp/webstart/	
	cp ../lib/substance.jar /tmp/webstart/
	cp -rf ../bin/* /tmp/webstart/
	
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
	cd /tmp/webstart; jar cf matheonF1webstart.jar * 

sign:
	cd /tmp/webstart; jarsigner -storepass hallo123 matheonF1webstart.jar sechel

access:
	cp -vf /tmp/webstart/matheonF1webstart.jar ~sechel/WWW/webstart/
	chmod og+r ~sechel/WWW/webstart/matheonF1webstart.jar
	cp -vf *.jnlp ~sechel/WWW/webstart/
	cp -vf *.png ~sechel/WWW/webstart/
	chmod og+r ~sechel/WWW/webstart/*.jnlp
	chmod og+r ~sechel/WWW/webstart/*.png
	rm -rf /tmp/webstart
		
