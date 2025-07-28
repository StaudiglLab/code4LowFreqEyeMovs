---------------
Installation and Dependencies
---------------

Dependencies:

MATLAB (R2020A) and the following toolboxes:
fieldtrip (20220617)
cbrewer2 (1.0.0.1)
daviolin (3.2.7)
dabarplot (3.2.4)
plot_areaerrorbar (1.3.1)

NOTE: 
1. Some functions may not be supported in MATLAB versions earlier than R2020a.
2. The complete analysis pipeline requires the full dataset located in \datafiles, which will be uploaded upon publication.
3. The results from the analyses are saved in \data4figs, which can be loaded now for figure generation.

---------------
Code to generate figures of the manuscript
---------------

The following code loads the saved outputs from the corresponding analyses and plots the figures

----------------------------------------- Figure 1 ---------------------------------------------
-- Fig1A: Beh\Beh_Exp1.m
-- Fig1B_left & Fig1D: Gaze\GazeAnalysis_Exp1.m
-- Fig1B_right: Saccades\SaccadesAnalysis_Exp1.m
-- Fig1C_left & Fig1E_left: TFR\TFR_Exp1.m
-- Fig1C_right & Fig1E_right: Source\Source_Exp1.m

----------------------------------------- Figure 2 ---------------------------------------------
-- Fig2A: Saccades\SaccadesAnalysis_Exp2.m
-- Fig2B: Source\Electrode_Location_Exp2.m
-- Fig2C_left & Fig2D_left: TFR\TFR_Exp2.m
-- Fig2C_right & Fig2D_right: Gaze\GazeAnalysis_Exp2.m
                      
----------------------------------------- Figure 3 ---------------------------------------------
-- Fig3A: Beh\Beh_Exp3.m
-- Fig3B_left & Fig3D_top: Saccades\SaccadesAnalysis_Exp3.m
-- Fig3B_right & Fig3D_bottom: Gaze\GazeAnalysis_Exp3.m
-- Fig3C_left & Fig3E_left: TFR\TFR_Exp3.m
-- Fig3E_right: Source\Source_Exp3.m

----------------------------------------- Figure 4 ---------------------------------------------
-- Fig4A: Beh\Beh_Exp4.m
-- Fig4B: Gaze\GazeAnalysis_Exp4.m
-- Fig4C: Saccades\SaccadesAnalysis_Exp4.m
-- Fig4D_left & Fig4E_left: TFR\TFR_Exp4.m
-- Fig4D_right & Fig4E_right: TFR\TFR_Exp4.m
-- Fig4D_top & Fig4E_top: Source\Source_Exp4.m

----------------------------------------- Figure 5 ---------------------------------------------
-- Fig5A: AlphaBeta\AlphaBeta_Saccade_Exp1.m
-- Fig5B: AlphaBeta\AlphaBeta_Saccade_Exp2.m
-- Fig5C: AlphaBeta\AlphaBeta_Saccade_Exp3.m
-- Fig5D: AlphaBeta\AlphaBeta_Saccade_Exp4.m

----------------------------------------- Figure 6 ---------------------------------------------
-- Fig6A: AlphaBeta\AlphaBeta_AfterSaccade_Exp1.m
-- Fig6B: AlphaBeta\AlphaBeta_AfterSaccade_Exp2.m
-- Fig6C: AlphaBeta\AlphaBeta_AfterSaccade_Exp3.m
-- Fig6D: AlphaBeta\AlphaBeta_AfterSaccade_Exp4.m

----------------------------------- Supplementary Figure 1 -------------------------------------
-- SuppleFig1: Saccades\Saccade_metrics.m

----------------------------------- Supplementary Figure 2 -------------------------------------
-- SuppleFig2: Gaze\ExplorationIndex_Example.m

----------------------------------- Supplementary Figure 3 -------------------------------------
-- SuppleFig3A: AlphaBeta\AlphaBeta_Saccade_Exp1.m
-- SuppleFig3B: AlphaBeta\AlphaBeta_Saccade_Exp2.m
-- SuppleFig3C: AlphaBeta\AlphaBeta_Saccade_Exp3.m
-- SuppleFig3D: AlphaBeta\AlphaBeta_Saccade_Exp4.m