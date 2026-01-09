function PrintTable(rowinf,colinf,tabinf)

    fprintf(rowinf{1},'');
    for j=1:size(colinf{2},2), fprintf(colinf{1},colinf{2}{j}); end;
    fprintf('\n');

    for r=1:size(rowinf{2},2)
        fprintf(rowinf{1},rowinf{2}{r});
        for j=1:size(colinf{2},2), fprintf(colinf{3}{j},tabinf(r,j)); end;
        fprintf('\n');
    end;
