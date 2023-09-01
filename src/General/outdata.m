mkdir('output')
filename = strsplit(dataloc1,{'\','/','.'});
dlmwrite('output/out_'+string(filename(end-1))+'.dat',opt_control_pts, 'delimiter',' ')
