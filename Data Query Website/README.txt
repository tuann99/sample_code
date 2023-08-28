* ABOUT *

Welcome to my final project!

This is a web-interface that allows users to search for anti-cancer drugs and get back information such as IC50, protein targets, and more.

The website is hosted at: http://bfx3.aap.jhu.edu/tnguy256/final/index.html

This tool utilizes data from the Wellcome Sanger Institute, specifically their project called, "Genomics of Drug Sensitivity in Cancer," or GDSC. Here is a link to their website: https://www.cancerrxgene.org/

I have created a SQL database (final.db), which will be queried using a CGI script, once the user submits their search term. I have tried my best to include both the trade name and generic name for as many drugs as I could, but this may not include the drug you're looking for.

In cases where your drug is not in the database, the webpage will still load, but the tables will not populate.

* TOOL USAGE *

This is a very simple tool, that only requires the user to type in a drug into the search bar, and either press the "submit" button, or press the "enter" key on the keyboard. Once submitted, the CGI script will run, and render the page with your results.

Example drugs to search:
Rapamycin
Bleomycin
Foretinib