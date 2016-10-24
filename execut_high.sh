#!/bin/tcsh -f

setenv dir `pwd`
#cd /apps/root/PRO/root/
cd /apps/root/5.34.05/root/
source bin/thisroot.csh
cd $dir

echo "2" > inp
echo "2.039" >> inp
echo "5" >> inp
@ z = $1 + 4
foreach k (`seq  $1 $z`)
echo "nt10_Aug16_newFrad_${k}.root" >> inp
end
echo "out.root" >> inp
./h10tot21<inp
