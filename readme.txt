CellectSeq:
===========
 In Silico Discovery of Antibodies Targeting Integral Membrane Proteins
 Combining In Situ Selections and Next-Generation Sequencing
 Published in Nature Communications Biology 2021 (in press)

Abstract:
=========
 Synthetic antibody (Ab) technologies are efficient and cost-effective platforms for the generation of
 monoclonal Abs against human antigens. Yet, they typically depend on purified proteins, which exclude
 integral membrane proteins that require the lipid bilayers to support their native structure and function.
 Here, we present an Ab discovery strategy, termed CellectSeq, for targeting integral membrane proteins
 on native cell surfaces in their complex environment. As proof of concept, we targeted three transmembrane
 proteins linked to cancer, tetraspanin CD151, carbonic anhydrase 9,and integrin-α11. First, we performed
 in situ cell-based selections against a 10^10 diversity library of synthetic phage displayed antigen-binding
 fragments (Fabs) (four CDRs were diversified: L3, H1, H2, and H3) to enrich synthetic Ab pools for antigen
 specific binders. Then, we designed next-generation sequencing procedures to explore Ab diversities and
 abundances in output selection pools. Finally, we developed motif-based scoring and sequencing error-filtering
 algorithms the comprehensive interrogation of next-generation sequencing pools to identify Abs with high
 diversities and specificities,  at extremely low abundances, which are impossible to identify using manual
 sampling or sequence abundances.

Strategy:
=========
 CellectSeq motif-based algorithm is based on the premise that  highly selective Abs are enriched with paratope
 motifs i.e., linear information) that recognize the target antigen, whereas non-selective Abs lack such enrichment.
 Therefore,  each Ab clone in the positive NGS pool we explore the entire space of linear information by exhaustively
 enumerating  possible motifs matching its CDR sequences, and obtain the frequencies of every motif in the positive
 and negative pools. According to the premise above, the high enrichment of  motifs in the positive pool relative to
 the negative pool implies the Ab candidate is potentially highly selective. , we analyzed each Ab in the positive pool
 for the selective binding to the target by scoring the separation between  two distributions of frequencies of the
 motifs in the positive and negative pools. To this end, we calculate the t-test  score the separation of the two
 distributions, then we calculate the p-value to evaluate the statistical significance  the t-test. Thus, the lower
 the p-value, the higher the separation between the two distributions, and consequently, the  the selectivity of the
 candidate Ab.
 
 CellectSeq works with NGS data of antibody sequences with up to four Complementarity-Determining Regions (CDRs),
 Also, it is able to handle multiple positive pools versus multiple negative pools, which make CellectSeq able to
 find highly selective antibodies with multi-specificities for multiple target antigens.
 

Folders:
========
 Source  : CellectSeq C++ source code
 INPUT   : Input files for NGS analysis
           "demux?.dat" defines NGS pools
           "tests?.dat" defines positive and negative pools
 RESULT  : NGS data and analysis for CD151, CA9, and ITGA11
 TEMP    : Temporary folder required during execution
 STEP01  : NGS equences pools
 STEP02  : Concatenated pools
 STEP03  : Motifs discovered in NGS pools
 STEP04  : Predicted high specificity sequences

Requirement:
============
 Operating system : Windows 10 (64 bits)
 Prerequisites    : Install cydwin from www.cygwin.com with default packages
                    Add cygwin and cygwin\bin directories to PATH environment variable
                    Download Boost C++ library from https://www.boost.org to CellectSeq main folder
                    Install Visual Studio C++ Comunity 2019 from https://visualstudio.microsoft.com
 Execution        : Open Visual Studio C++ Community 2019
                    Browse to CellectSeq folder
                    Open CellectSeq project
                    Compile then execute

Development:
============
Abdellali Kelil, PhD
Reseach Associate and Senior Bioinformatician
Sidhu Laboratory - Donnelly Centre - University of Toronto
160 College Street Room 810 - Toronto, Ontario CANADA M5S 3E1
abdellali.kelil@utoronto.ca

Acknowledgment:
===============
 I thank Eugenio Gallo and Jarrett Adams from Toronto Recombinant Antibody Centre for their help in developping SelelctSeq.
 I thank Dax Torti and Tanja Durbic from the Donnelly Sequencing Centre for their help with NGS runs.
 I thank Meghan McLaughlin and Wei Ye for designing the primers used in NGS reads and improvements of cellular selections.
 I thank Sachdev Sidhu from the Donnelly Centre for supervising and financing the whole project.
