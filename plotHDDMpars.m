function plotHDDMpars()
% csv = readtable('/Users/kloosterman/Dropbox/PROJECTS/EyeMem/HDDM/EyeMem_hddm_test.csv');
csv = readtable('/Users/kloosterman/Dropbox/PROJECTS/COBRA/hddm/123back_bias/params_HDDMbias.csv')
csv(1:10,:)

ind = find(strcmp(csv.Var1, 'z(1-back)'));
csv(ind:ind+2,:)
ind = find(strcmp(csv.Var1, 'dc(1-back)'));
csv(ind:ind+2,:)
ind = find(strcmp(csv.Var1, 'v(1-back)'));
csv(ind:ind+2,:)
ind = find(strcmp(csv.Var1, 'a(1-back)'));
csv(ind:ind+2,:)
ind = find(strcmp(csv.Var1, 't(1-back)'));
csv(ind:ind+2,:)

figure