#!/bin/bash
folderin="./FT/"
folderout="./"

folderimage="../Images_231018/"
force=(5 6 7 8 9 10 11 12)
Dt=(10 20 40 80 160 320)

for f in ${force[*]}
do
	fileout="${folderout}${f}_FT.gpl"
	echo "
unset multiplot
clear
reset
set print 'A_force.txt'
set key bot right
set term qt
set xra[0:16]
set yra[0:7]
set xla '{/Symbol D}x/<{/Symbol D}x>'
set yla 'log(p({/Symbol D}x)/p(-{/Symbol D}x))/<{/Symbol D}x>'
set title 'f=$f pN'
"> $fileout
	for ((k=0;$k<6;k++))
	do
		file="${folderin}${f}_${Dt[$k]}FT.txt"
		if [[ "$k" != "0" ]]
		then
		echo "rep '$file' u 1:2 ps 2 lw 3 t '{/Symbol D}t=${Dt[$k]}dt'">>${fileout}
		fi
		if [[ "$k" == "0" ]]
		then
		echo "p '$file' u 1:2 ps 2 lw 3 t '{/Symbol D}t=${Dt[$k]}dt'">>${fileout}
		fi
	done
	file="${folderin}${f}_${Dt[$4]}FT.txt"
	echo "f(x)=A*x+B
B=0.0
fit f(x) '${file}' via A
fitit=sprintf('A=%.2f', A)
print ${f},A
unset print
rep f(x) lw 3 t fitit
set term post eps enhanced color solid 24
set output '${folderimage}${f}_FT.eps'
rep
unset output
set term qt">>${fileout}
done
