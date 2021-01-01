
# Code Policy
Please state “Cell Flux Model” and “Keisuke Inomura” in the acknowledgement when your
publication includes the results based on the original/revised code. Or you may consider
including Keisuke Inomura as a co-author depending on the contribution. In either case, the
publication must cite the following papers:
Masuda T, Inomura K, Takahata N, Shiozaki T, Sano Y, Deutsch C, Prášil O, Furuya K. 2020. Heterogeneous nitrogen fixation rates confer energetic advantage and expanded ecological niche of unicellular diazotroph populations. Communications Biology 3:172.
Inomura K, Bragg J, Riemann L, Follows MJ. 2018. A quantitative model of nitrogen fixaion in the presence of ammonium. PLoS ONE 13:e0208282.
Inomura K, Bragg J, Follows MJ. 2017. A quantitative analysis of the direct and indirect costs of nitrogen fixation: a model based on Azotobacter vinelandii. ISME Journal 11:166–175.
(Papers downloaded from https://www.inomura.com/papers)

Keisuke Inomura (University of Rhode Island)
inomura@uri.edu


# CrocoHetero
Simulating Heterogenaity of Crocosphaera
**************
Preparation:
Download all the files.
In each file, path must be adjusted to your own folders.
**************
Fig. 6a,b
Run 617_04_20_09_ColorRevision.py

Fig. 6c,d
1. Run 617_04_22_06_Stack_Dd_normalized_by_Vch_T26.py with "CreatingNormal = 0", which will produce the result with Rv = 1
2. Run the same file with "CreatingNormal = 1"
* make sure to adjust the path.

Fig. 7a,b
1. Run 631_00_06_2D_Rn_Ef_HighRes02.py
2. For plotting, run 641_00_09_2D_optRn_N2FixMax.py
* make sure to adjust the path.

Fig. 7c,d
1. Run 617_04_22_04_Stack_Dd_normalized_by_Vch_T26.py with "CreatingNormal = 0", which will produce the result with Rv = 1
2. Run 633_00_07_2D_Dd_Rn_HighRes02_T26.py
3. Run 643_00_11_OptRn.py
* make sure to adjust the path.

Fig. 8
Run L001_01_04_read_Matlab_file.py
* make sure to adjust the path.

Fig. S3
Run 617_04_20_09_ColorRevision.py with "Ef=0.1"

Fig. S4
Run L003_00_01_Shiozaki2018data.py

Fig. S5
Run L002_00_00_EN0.2_0.5_compare.py
