Instructions on using stellar cluster distance estimate pipeline
  - Code by Mitchel Miners, Ali Bo Khamseen, and Maximillian Mueth

--OVERVIEW--
The pipeline reads in queries of GAIA and SIMBAD data (w/ optional cluster membership list) and returns distance estimates from trigonometric parallax, spectroscopic parallax, main sequence (MS) fitting, variable extinction analysis, and the cepheid PL relation (untested).
Sample files for GAIA and SIMBAD query results are included for reference. The header of the SIMBAD results was altered slightly to make data parsing simplier. A sample GAIA query is included for reference on what data is used in our pipeline. SIMBAD queries are done using the database website's cone search.
Users set the query result filenames (and optional cluster membership filename) at the top of the pipeline.py file. A magnitude filter can also be adjusted and the total-to-selective extinction ratio R can be set to 3.1, the galactic average, instead of relying on calculated values in your cluster. The TestData option tweaks database crossmatching to ensure better cluster membership at exetremely close distances and generates a MS from a cubic fit to your cluster. The MS used in the zams() function was created using the Hyades cluster in TestData mode.

--USAGE STEPS--
1) Run queries in SIMBAD and GAIA for your cluster. Sample GAIA queries and results are included in this project path. SIMBAD queries can be done using the official website's cone search tool. Sample SIMBAD query results are included for reference. Note that the SIMBAD result header was manually modified for ease of access in python.
2) (OPTIONAL) Run membership determinations on your cluster using UPMASK. A sample file with membership probabilies in included for reference. Members with probabilites >=.75 are kept. If you do not wish to use a membership list, set didMembership to false in the user parameter section of pipeline.py and proper motions will be used to exclude members outside a specified number of standard deviations.
3) Set user parameters at the top of the pipeline.py file.
4) Place pipeline.py and all data files in the same file directory. Then run pipeline.py through the command line.
5) Several plots about your cluster should display along with distance and other information printed in the command line.
