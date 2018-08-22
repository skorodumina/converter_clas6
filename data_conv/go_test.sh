#!/bin/tcsh -f





@ x = 0
@ y = 1
while ($x <= 1989)
echo '#\!/bin/tcsh -f' > run_data$y.sh
echo ' ' >> run_data$y.sh

echo 'setenv dir `pwd`' >> run_data$y.sh
echo 'cd /apps/root/5.34.21/root/' >> run_data$y.sh
echo 'source bin/thisroot.csh' >> run_data$y.sh
echo 'cd $dir' >> run_data$y.sh
echo ' ' >> run_data$y.sh


@ z = $x + 1



sed -e "1,$z d" file_list |  sed  -e "7q" | sed -e 's$/mss$jget /mss$g' | sed -e "s/.root/.root ./g">> run_data$y.sh
echo ' ' >> run_data$y.sh


echo 'echo "1" > inp' >> run_data$y.sh
echo 'echo "2.039" >> inp' >> run_data$y.sh
echo 'echo "7" >> inp' >> run_data$y.sh
sed -e "1,$z d" file_list |  sed  -e "7q" | sed -e 's$/mss/clas/e1e/production/pass1/h10/$$g' | sed -e 's$a1ntp$echo "a1ntp$' | sed -e 's$root$root" >> inp$'>> run_data$y.sh
echo 'echo "out.root" >> inp' >> run_data$y.sh
echo './h10tot21<inp' >> run_data$y.sh
chmod +x run_data$y.sh
@ y += 1

@ x += 7
end
echo $y

@ s = $y - 2
echo $s
