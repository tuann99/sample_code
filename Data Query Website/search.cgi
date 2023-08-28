#!/usr/local/bin/python3

import cgi
import jinja2 as j
import sqlite3 as sql

# get the search term, connect to db, load template
form = cgi.FieldStorage()
get = form.getvalue('drug')
# want all input to be uniform so lower/upper case doesn't break code
search_term = str(get).lower()  
drug_name = str(search_term)

templateLoader = j.FileSystemLoader(searchpath="/var/www/html/tnguy256/final/html")
env = j.Environment(loader=templateLoader)
template = env.get_template('template.html')

conn = sql.connect('/var/www/html/tnguy256/final/final.db')
curs = conn.cursor()

# query db
query = """
SELECT t.drug_target, p.drug_target_pathway, i.ic50, c.cell_line_name
FROM target t
JOIN keys k ON t.target_id = k.target_id
JOIN pathway p ON k.path_id = p.path_id
JOIN ic50 i ON k.entry_id = i.entry_id
JOIN drug_name_alt a ON k.drug_id = a.drug_id
JOIN drug d ON k.drug_id = d.drug_id
JOIN cell_line c ON k.cell_line_id = c.cell_line_id
WHERE d.drug_name LIKE ? OR a.drug_name_alt LIKE ?
GROUP BY t.drug_target, p.drug_target_pathway, c.cell_line_name;
"""
data = curs.execute(query, (search_term, search_term))
rows = data.fetchall()

# make list to append results to
entries = list()

for row in rows: 
    entries.append(
        {'drug_target':row[0],
         'drug_target_pathway':row[1],
         'ic50':row[2],
         'cell_line':row[3]
         })

# render html
print("Content-Type: text/html\n\n")
print(template.render(sql_query=entries, drug_name=drug_name))

curs.close()
conn.close()