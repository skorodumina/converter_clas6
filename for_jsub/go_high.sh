#!/bin/tcsh -f
#setenv size 0
foreach k (`seq  1 5 75000`)


@ z = $k / 5 + 1
#echo $z
#setenv size `ls -all  /cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${z}.root | grep size | sed -e 's/size=//g'`
#setenv size `ls -l $/cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${z}.root | awk '{print $5}'| sed -e 's/size=//g'`
#333333setenv size `stat -c %s $/cache/mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${z}.root`


echo "PROJECT: e1e" > jsub_new
echo "TRACK: simulation" >>jsub_new
echo "JOBNAME: convert_skor_fermi_${k}" >>jsub_new
echo "MAIL: skorodum@jlab.org" >>jsub_new
echo "SINGLE_JOB: TRUE" >>jsub_new
echo "MEMORY: 3000 MB" >>jsub_new
echo "DISK_SPACE: 15 GB" >>jsub_new
echo "OS:  centos65" >>jsub_new
@ k1 = $k + 1
@ k2 = $k + 2
@ k3 = $k + 3
@ k4 = $k + 4
echo "INPUT_FILES: /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/nt10/nt10_Aug16_newFrad_${k}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/nt10/nt10_Aug16_newFrad_${k1}.root  /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/nt10/nt10_Aug16_newFrad_${k2}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/nt10/nt10_Aug16_newFrad_${k3}.root /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/nt10/nt10_Aug16_newFrad_${k4}.root /volatile/clas/clase1-6/skorodum/converter_gleb/h10tot21" >>jsub_new
echo "COMMAND: /volatile/clas/clase1-6/skorodum/converter_gleb/big_sim_conv_25May2016/execut_high.sh ${k}" >>jsub_new
echo "OUTPUT_DATA: out.root" >>jsub_new
echo "OUTPUT_TEMPLATE: /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${z}.root" >>jsub_new

#if ($size < 220000000) then
if (!(-e /mss/clas/e1e/production/simulation_2pi/sim_skorodum_Aug2016/converted_Aug2016/out_newFrad_Aug2016_${z}.root )) then
/site/bin/jsub jsub_new
echo ${z}
endif

rm jsub_new

end
