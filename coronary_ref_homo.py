import sqlite3
from sqlite3 import Error
from pathlib import Path


class CoronaryRefHomo:

    rsid_map:dict[str, dict] = {}


    def setup(self, parent, result_cursor:sqlite3.Cursor, data_cursor:sqlite3.Cursor, sql_insert:str) -> None:
        self.parent = parent
        self.sql_insert = sql_insert
        self.result_cursor: sqlite3.Cursor = result_cursor
        self.data_cursor: sqlite3.Cursor = data_cursor

        sql:str = "SELECT rsID, Ref_allele FROM coronary_disease WHERE state = 'ref' AND zygosity = 'hom'"
        self.data_cursor.execute(sql)
        rows:list = self.data_cursor.fetchall()
        for rsid, ref_allele in rows:
            self.rsid_map[rsid] = {'exist':True, 'ref':ref_allele}


    def process_row(self, row):
        rsid:str = str(row['dbsnp__rsid'])
        if rsid == '':
            return

        if not rsid.startswith('rs'):
            rsid = "rs"+rsid

        item:dict = self.rsid_map.get(rsid)
        if item:
            self.rsid_map[rsid]['exist'] = False


    def end(self) -> None:
        for rsid in self.rsid_map:
            if self.rsid_map[rsid]['exist']:
                ref:str = self.rsid_map[rsid]['ref']
                ref = ref+ref

                query:str = "SELECT Risk_allele, Gene, Genotype, Conclusion, Weight, PMID, Population, GWAS_study_design, P_value " \
                        f"FROM coronary_disease WHERE rsID = '{rsid}' AND Genotype = '{ref}';"

                self.data_cursor.execute(query)
                row:list = self.data_cursor.fetchone()
                if row:
                    task:tuple = (rsid, row[1], row[0], row[0]+"/"+row[0], row[3], row[4], row[5], row[6], row[7], row[8],
                            self.parent.get_color(row[4], 0.6))
                    self.result_cursor.execute(self.sql_insert, task)