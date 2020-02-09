function data = choseData(data)
    if data == 1
        disp('import data: gse2403_BCLL_norm.csv')
        file = importdata('../data/gse2403_BCLL_norm.csv');
        data = file.data; %size 18125x21
    
    elseif data == 2
        disp('import data: gds2350_YC_norm.csv')
        file = importdata('../data/gds2350_YC_norm.csv');
        data = file.data; %size 5847x50
                
    elseif data == 3
        disp('import data: gse952_RAT_norm.csv')
        file = importdata('../data/gse952_RAT_norm.csv');
        data = file.data; %size 7751x122
        
    elseif data == 4
        disp('import data: gse2034_PBC_norm.csv')
        file = importdata('../data/gse2034_PBC_norm.csv');
        data = file.data; %size 21225x286
    
    end
end