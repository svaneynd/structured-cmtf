# structured-cmtf
Pipeline for structured coupled matrix-tensor factorization (CMTF) of simultaneously recorded EEG and fMRI data.

This directory contains a collection of all MATLAB functions that are needed to perform the computations described in

Van Eyndhoven S., Dupont P., Tousseyn S., Vervliet N., Van Paesschen W., Van Huffel S., Hunyadi B., "Augmenting interictal mapping with neurovascular coupling biomarkers by structured factorization of epileptic EEG and fMRI data", Internal Report 20-59, ESAT-STADIUS, KU Leuven (Leuven, Belgium), 2020., Lirias number: x.

If you use (parts of) this code, please cite this work.

To get started: execute 'main_pipeline.m'. Providing your own data is necessary. However, the code can generate toy data, which allow to run the pipeline end-to-end, and observe the type of outputs that are being generated.

In order to run the code, the user needs 1) an installation of the [Statistical Parametric Mapping (SPM) toolbox](https://www.fil.ion.ucl.ac.uk/spm/software/), version 8 or 12; 2) the mwf-artifact-removal toolbox:

[1] Somers, B., Francart, T. and Bertrand, A. (2018). A generic EEG artifact removal algorithm based on the multi-channel Wiener filter. Journal of Neural Engineering, 15(3), 036007. DOI: 10.1088/1741-2552/aaac92 

[2] Somers, B., Francart, T. and Bertrand, A. (2017). MWF toolbox for EEG artifact removal. Available online, URL: http://www.github.com/exporl/mwf-artifact-removal


For convenience, we include partial or complete code of less common toolboxes, which our code relies on:

Tensorlab: 

[3] Vervliet N., Debals O., Sorber L., Van Barel M. and De Lathauwer L. Tensorlab 3.0, Available online, Mar. 2016. URL: https://www.tensorlab.net/

Brainnetome:

[4] Fan, L., Li, H., Zhuo, J., Zhang, Y., Wang, J., Chen, L., Yang, Z., Chu, C., Xie, S., Laird, A.R., Fox, P.T., Eickhoff, S.B., Yu, C. & Jiang, T. The Human Brainnetome Atlas: A New Brain Atlas Based on Connectional Architecture. Cerebral Cortex, 26 (8): 3508-3526,(2016). 

[5] URL: http://atlas.brainnetome.org

CONN:

[6] Whitfield-Gabrieli, S., & Nieto-Castanon, A. (2012). Conn: A functional connectivity toolbox for correlated and anticorrelated brain networks. Brain connectivity, 2(3), 125-141

[7] CONN toolbox (www.nitrc.org/projects/conn, RRID:SCR_009550)

kde

[8] Kernel density estimation via diffusion Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010) Annals of Statistics, Volume 38, Number 5, pages 2916-2957.

[9] URL: https://nl.mathworks.com/matlabcentral/fileexchange/58312-kernel-density-estimator-for-high-dimensions?s_tid=prof_contriblnk

Author: Simon Van Eyndhoven (simon.vaneyndhoven@gmail.com)
