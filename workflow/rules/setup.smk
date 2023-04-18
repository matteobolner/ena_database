rule split_db_by_unit:
    input:
        #db="data/suina/pafr_test.csv"
        db="data/suina/test.tsv"
    output:
        outfile="data/split_db/{sample}/{unit}.tsv"
    params:
        sample= lambda wc:wc.sample,
        unit= lambda wc:wc.unit
    run:
        tempdf=pd.read_table(input.db)
        tempdf=tempdf[(tempdf['accession']==params.sample)&(tempdf['run_accession']==params.unit)]
        tempdf.to_csv(output.outfile, index=False, sep='\t')
