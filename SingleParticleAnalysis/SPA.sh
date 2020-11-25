#!/bin/bash
echo 'enter start sample number:'
read nstart
echo 'enter stop sample number:'
read nstop
foldername='SingleParticleAnalysis'
exc='SingleParticleAnalysis.exe'
for ((i = $nstart; i <= $nstop; i++)); do
	echo $i' start'
	cd ..
	cp -r ./$foldername ./$i/
	cd ./$i/$foldername
	./$exc >'SPA.log' 2>&1
	cd ..
	echo $i' finshed'
done
