#!/bin/bash
export JAVA_HOME=/usr/site-local/jdk1.8
ant release
java -Xms4G -Xmx4G -classpath release/$(date +%F)/DiscreteConformalLab.jar de.varylab.discreteconformal.holomorphicformsexperiments.convergence.FibonacciClusteringG2 >> FibonacciClusteringG2.log 2>&1

