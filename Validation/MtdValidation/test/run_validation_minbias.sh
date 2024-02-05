cmsRun mtdValidation_cfg.py -minbias &> validation_minbias.log;
cmsRun mtdHarvesting_cfg.py -minbias &>> validation_minbias.log;
mv DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root DQM_MinBias.root