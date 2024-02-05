cmsRun mtdValidation_cfg.py &> validation_ttbar.log;
cmsRun mtdHarvesting_cfg.py &>> validation_ttbar.log;
mv DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root DQM_TTbar.root