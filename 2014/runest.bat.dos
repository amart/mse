rm pk13_3.par
rm pk13_3.rep

cp -p pk13_proj.dat pk13_proj.$1.$2.dat

./pk13_3 -nox -nohess -ind pk13_proj.dat

cp -p pk13_3.rep pk13_3.$1.$2.rep
cp -p pk13_3.par pk13_3.$1.$2.par

