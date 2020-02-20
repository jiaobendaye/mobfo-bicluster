rootPath = pwd;
datas = {'BCLL', 'PBC', 'RAT', 'YC'};
data_N_map = containers.Map;
data_N_map('BCLL') = 12185;
data_N_map('PBC') = 21225;
data_N_map('RAT') = 7751;
data_N_map('YC') = 5847;


for i=1:length(datas)
    data = char(datas(i));
    N = data_N_map(data);
    dataPath = char(fullfile(rootPath, data));
    dirOut = dir(char([dataPath, '/*.mat']));
    matFiles = {dirOut.name};
    matCnt = fix(length(matFiles) / 2);
    if matCnt == 0
        continue;
    end
    txtFile = [dataPath, '.txt'];
    disp(txtFile)
    fid = fopen(txtFile, 'w');
    for k=1:matCnt
        load(fullfile(dataPath, [num2str(k*20, '%04d'), '_', 'bic.mat']));
        load(fullfile(dataPath, [num2str(k*20, '%04d'), '_', 'scores.mat']));
        for l=1:20
            bic = pareto_set(l,:);
            scores = Pscores(l,:);
            [~, genes] = find(bic(1:N)==1);
            [~, conditions] = find(bic(N+1:length(bic))==1);
            %matlab index is start from 1, but python is 0
            genes = genes -1;
            conditions = conditions -1;

            fprintf(fid, '%f\n',scores(1));
            fprintf(fid, '%d\n',scores(3));
            fprintf(fid, '%d\n',scores(4));
            fprintf(fid, '%f\n',scores(5));
            fprintf(fid, '%s\n',num2str(genes));
            fprintf(fid, '%s\n',num2str(conditions));
        end
    end
    fclose(fid);
end