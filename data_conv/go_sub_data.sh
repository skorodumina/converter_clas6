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

foreach k (`seq  1 $s`)

echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: convert_skor_data_${k}" >>jsub_new
echo "MAIL: skorodum@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3000 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "OS:  centos7" >>jsub_new
echo "INPUT_FILES:/volatile/clas/clase1-6/skorodum/converter_gleb/data_conv/h10tot21" >>jsub_new
echo "COMMAND: /volatile/clas/clase1-6/skorodum/converter_gleb/data_conv/run_data${k}.sh" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/home/skorodum/e1e/data_2pi_conv_8Mar2018/out_conv_2pi_8Mar2018_${k}.root" >>jsub_new



if (!(-e /mss/home/skorodum/e1e/data_2pi_conv_8Mar2018/out_conv_2pi_8Mar2018_${k}.root)) then
echo $k
/site/bin/jsub jsub_new
endif

rm jsub_new

end
