from oakvar import BasePostAggregator
from pathlib import Path
import sys
cur_path = str(Path(__file__).parent)
sys.path.append(cur_path)
import sqlite3
import coronary_ref_homo


class CravatPostAggregator (BasePostAggregator):
    sql_insert:str = """ INSERT INTO coronary (
                        rsid,
                        gene,
                        risk,
                        genotype,
                        conclusion,
                        weight,
                        pmid,
                        population,
                        studydesign,
                        pvalue,
                        weightcolor        
                    ) VALUES (?,?,?,?,?,?,?,?,?,?,?) """
    ref_homo:coronary_ref_homo.CoronaryRefHomo = coronary_ref_homo.CoronaryRefHomo()

    def check(self):
        return True

    def setup (self):
        modules_path:str = str(Path(__file__).parent)
        sql_file:str = modules_path + "/data/coronary.sqlite"
        if Path(sql_file).exists():
            self.data_conn:sqlite3.Connection = sqlite3.connect(sql_file)
            self.data_cursor:sqlite3.Cursor = self.data_conn.cursor()

        self.result_path:Path = Path(self.output_dir, self.run_name + "_longevity.sqlite")
        self.result_conn:sqlite3.Connection = sqlite3.connect(self.result_path)
        self.result_cursor:sqlite3.Cursor = self.result_conn.cursor()
        sql_create:str = """ CREATE TABLE IF NOT EXISTS coronary (
            id integer NOT NULL PRIMARY KEY,
            rsid text,
            gene text,
            risk text,
            genotype text,
            conclusion text,
            weight float,
            pmid text,
            population text,
            studydesign text,
            pvalue text,
            weightcolor text            
            )"""
        self.result_cursor.execute(sql_create)
        self.result_conn.commit()
        self.result_cursor.execute("DELETE FROM coronary;")
        self.ref_homo.setup(self, self.result_cursor, self.data_cursor, self.sql_insert)

    
    def cleanup (self):
        if self.result_cursor is not None:
            self.result_cursor.close()
        if self.result_conn is not None:
            self.result_conn.commit()
            self.result_conn.close()
        if self.data_cursor is not None:
            self.data_cursor.close()
        if self.data_conn is not None:
            self.data_conn.close()
        return


    def get_color(self, w:float, scale:float = 1.5) -> str:
        w = float(w)
        if w < 0:
            w = w * -1
            w = 1 - w * scale
            if w < 0:
                w = 0
            color:str = format(int(w * 255), 'x')
            if len(color) == 1:
                color = "0" + color
            color = "ff" + color + color
        else:
            w = 1 - w * scale
            if w < 0:
                w = 0
            color = format(int(w * 255), 'x')
            if len(color) == 1:
                color = "0" + color
            color = color + "ff" + color

        return color

        
    def annotate (self, input_data):
        rsid:str = str(input_data['dbsnp__rsid'])
        if rsid == '':
            return

        self.ref_homo.process_row(input_data)

        if not rsid.startswith('rs'):
            rsid = "rs" + rsid

        alt:str = input_data['base__alt_base']
        ref:str = input_data['base__ref_base']

        query:str = "SELECT Risk_allele, Gene, Genotype, Conclusion, Weight, PMID, Population, GWAS_study_design, P_value " \
                f"FROM coronary_disease WHERE rsID = '{rsid}';"

        self.data_cursor.execute(query)
        rows:list[Any] = self.data_cursor.fetchall()

        if len(rows) == 0:
            return

        zygot:str = input_data['vcfinfo__zygosity']
        genome:str = alt + "/" + ref
        gen_set:set[str] = {alt, ref}
        if zygot == 'hom':
            genome = alt + "/" + alt
            gen_set = {alt, alt}
        for row in rows:
            allele = row[0]
            row_gen:set[str] = {row[2][0], row[2][1]}

            if gen_set == row_gen:
                task:tuple = (rsid, row[1], allele, genome, row[3], float(row[4]), row[5], row[6], row[7],
                        row[8], self.get_color(row[4], 0.6))
                self.result_cursor.execute(self.sql_insert, task)

        return {"col1":""}


    def postprocess(self):
        self.ref_homo.end()
        pass
