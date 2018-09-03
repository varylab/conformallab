#!/bin/bash
export JAVA_HOME=/usr/site-local/jdk1.8
ant release
java -Xms4G -Xmx4G -classpath release/$(date +%F)/DiscreteConformalLab.jar de.varylab.discreteconformal.holomorphicformsexperiments.convergence.LawsonG2Scattering >> LawsonG2-Scatter.log 2>&1

