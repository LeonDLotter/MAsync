% get go descriptions from gofigure results via ABAnnotate

gofigure = readtable('biological_process_GCEA_ale_z_GO.tsv', 'FileType', 'delimitedtext');

n_terms = 5;

descs = {};
for i=1:n_terms
    terms = gofigure.members{i};
    terms = erase(terms, ["['", "']"]);
    terms = split(terms, "', '");
    descs{i,1} = get_go_desc(terms);
end

descs_str = cellfun(@(x) strjoin(string(x), ', '), descs);
writetable(array2table(descs_str), 'gofigure_res_descriptions_5.csv');
