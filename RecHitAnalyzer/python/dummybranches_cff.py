import FWCore.ParameterSet.Config as cms
#tableName = cms.string('/frozen/2018/2e34/v3.2/HLT/V1')

dummy_branches_for_PbPb_2018_HLT = cms.untracked.vstring([
    'DST_Physics_v7',
    'DST_ZeroBias_v2',
    'HLT_AK8PFHT750_TrimMass50_v12',
    'HLT_AK8PFHT800_TrimMass50_v12',
    'HLT_AK8PFHT850_TrimMass50_v11',
    'HLT_AK8PFHT900_TrimMass50_v11',
    'HLT_AK8PFJet140_v15',
    'HLT_AK8PFJet15_v3',
    'HLT_AK8PFJet200_v15',
    'HLT_AK8PFJet25_v3',
    'HLT_AK8PFJet260_v16',
    'HLT_AK8PFJet320_v16',
    'HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p17_v2',
    'HLT_AK8PFJet330_TrimMass30_PFAK8BTagDeepCSV_p1_v2',
    'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np2_v2',
    'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_v2',
    'HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_p02_v3',
    'HLT_AK8PFJet360_TrimMass30_v18',
    'HLT_AK8PFJet380_TrimMass30_v11',
    'HLT_AK8PFJet400_TrimMass30_v12',
    'HLT_AK8PFJet400_v16',
    'HLT_AK8PFJet40_v16',
    'HLT_AK8PFJet420_TrimMass30_v11',
    'HLT_AK8PFJet450_v16',
    'HLT_AK8PFJet500_v16',
    'HLT_AK8PFJet550_v11',
    'HLT_AK8PFJet60_v15',
    'HLT_AK8PFJet80_v15',
    'HLT_AK8PFJetFwd140_v14',
    'HLT_AK8PFJetFwd15_v3',
    'HLT_AK8PFJetFwd200_v14',
    'HLT_AK8PFJetFwd25_v3',
    'HLT_AK8PFJetFwd260_v15',
    'HLT_AK8PFJetFwd320_v15',
    'HLT_AK8PFJetFwd400_v15',
    'HLT_AK8PFJetFwd40_v15',
    'HLT_AK8PFJetFwd450_v15',
    'HLT_AK8PFJetFwd500_v15',
    'HLT_AK8PFJetFwd60_v14',
    'HLT_AK8PFJetFwd80_v14',
    'HLT_CaloJet500_NoJetID_v12',
    'HLT_CaloJet550_NoJetID_v7',
    'HLT_DiPFJetAve100_HFJEC_v16',
    'HLT_DiPFJetAve140_v13',
    'HLT_DiPFJetAve160_HFJEC_v16',
    'HLT_DiPFJetAve200_v13',
    'HLT_DiPFJetAve220_HFJEC_v16',
    'HLT_DiPFJetAve260_v14',
    'HLT_DiPFJetAve300_HFJEC_v16',
    'HLT_DiPFJetAve320_v14',
    'HLT_DiPFJetAve400_v14',
    'HLT_DiPFJetAve40_v14',
    'HLT_DiPFJetAve500_v14',
    'HLT_DiPFJetAve60_HFJEC_v15',
    'HLT_DiPFJetAve60_v14',
    'HLT_DiPFJetAve80_HFJEC_v16',
    'HLT_DiPFJetAve80_v13',
    'HLT_DoublePFJets100_CaloBTagDeepCSV_p71_v2',
    'HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v2',
    'HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v2',
    'HLT_DoublePFJets200_CaloBTagDeepCSV_p71_v2',
    'HLT_DoublePFJets350_CaloBTagDeepCSV_p71_v2',
    'HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v2',
    'HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_v2',
    'HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_v2',
    'HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_v2',
    'HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v2',
    'HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v2',
    'HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v2',
    'HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v2',
    'HLT_PFHT1050_v18',
    'HLT_PFHT180_v17',
    'HLT_PFHT250_v17',
    'HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v3',
    'HLT_PFHT330PT30_QuadPFJet_75_60_45_40_v9',
    'HLT_PFHT350MinPFJet15_v9',
    'HLT_PFHT350_v19',
    'HLT_PFHT370_v17',
    'HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_v8',
    'HLT_PFHT400_SixPFJet32_v8',
    'HLT_PFHT430_v17',
    'HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_v7',
    'HLT_PFHT450_SixPFJet36_v7',
    'HLT_PFHT500_PFMET100_PFMHT100_IDTight_v12',
    'HLT_PFHT500_PFMET110_PFMHT110_IDTight_v12',
    'HLT_PFHT510_v17',
    'HLT_PFHT590_v17',
    'HLT_PFHT680_v17',
    'HLT_PFHT700_PFMET85_PFMHT85_IDTight_v12',
    'HLT_PFHT700_PFMET95_PFMHT95_IDTight_v12',
    'HLT_PFHT780_v17',
    'HLT_PFHT800_PFMET75_PFMHT75_IDTight_v12',
    'HLT_PFHT800_PFMET85_PFMHT85_IDTight_v12',
    'HLT_PFHT890_v17',
    'HLT_PFJet140_v19',
    'HLT_PFJet15_v3',
    'HLT_PFJet200_v19',
    'HLT_PFJet25_v3',
    'HLT_PFJet260_v20',
    'HLT_PFJet320_v20',
    'HLT_PFJet400_v20',
    'HLT_PFJet40_v21',
    'HLT_PFJet450_v21',
    'HLT_PFJet500_v21',
    'HLT_PFJet550_v11',
    'HLT_PFJet60_v21',
    'HLT_PFJet80_v20',
    'HLT_PFJetFwd140_v18',
    'HLT_PFJetFwd15_v3',
    'HLT_PFJetFwd200_v18',
    'HLT_PFJetFwd25_v3',
    'HLT_PFJetFwd260_v19',
    'HLT_PFJetFwd320_v19',
    'HLT_PFJetFwd400_v19',
    'HLT_PFJetFwd40_v19',
    'HLT_PFJetFwd450_v19',
    'HLT_PFJetFwd500_v19',
    'HLT_PFJetFwd60_v19',
    'HLT_PFJetFwd80_v18',
    'HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v8',
    'HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v8',
    'HLT_QuadPFJet103_88_75_15_v5',
    'HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v8',
    'HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v8',
    'HLT_QuadPFJet105_88_76_15_v5',
    'HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v8',
    'HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v8',
    'HLT_QuadPFJet111_90_80_15_v5',
    'HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v8',
    'HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_v8',
    'HLT_QuadPFJet98_83_71_15_v5',
    'HLT_Rsq0p35_v15',
    'HLT_Rsq0p40_v15',
    'HLT_RsqMR300_Rsq0p09_MR200_4jet_v15',
    'HLT_RsqMR300_Rsq0p09_MR200_v15',
    'HLT_RsqMR320_Rsq0p09_MR200_4jet_v15',
    'HLT_RsqMR320_Rsq0p09_MR200_v15',
    'HLT_SingleJet30_Mu12_SinglePFJet40_v11',
    'HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_v1',
    'HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4',
    'HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_v1',
    'HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v1',
    'HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_v1',
    'HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v1',
    'HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_v1',
    'HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v1',
    'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v11',
    'HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12',
    'HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_v12',
    'HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_v12',
    'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12',
    'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v8',
    'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_v8',
    'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v8',
    'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140_v3',
    'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v12',
    'HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v12',
    'HLT_Photon35_TwoProngs35_v1',
    'HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_v1',
    'HLT_VBF_DoubleMediumChargedIsoPFTauHPS20_Trk1_eta2p1_v1',
    'HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_v1'
    ])

dummy_branches_for_PbPb_2018_L1 = cms.untracked.vstring([
    'L1_AlwaysTrue',
    'L1_BPTX_AND_Ref1_VME',
    'L1_BPTX_AND_Ref3_VME',
    'L1_BPTX_AND_Ref4_VME',
    'L1_BPTX_BeamGas_B1_VME',
    'L1_BPTX_BeamGas_B2_VME',
    'L1_BPTX_BeamGas_Ref1_VME',
    'L1_BPTX_BeamGas_Ref2_VME',
    'L1_BPTX_NotOR_VME',
    'L1_BPTX_OR_Ref3_VME',
    'L1_BPTX_OR_Ref4_VME',
    'L1_BPTX_RefAND_VME',
    'L1_BptxMinus',
    'L1_BptxMinus_NotBptxPlus',
    'L1_BptxOR',
    'L1_BptxPlus',
    'L1_BptxPlus_NotBptxMinus',
    'L1_BptxXOR',
    'L1_Castor1',
    'L1_CastorHighJet',
    'L1_CastorHighJet_BptxAND',
    'L1_CastorHighJet_MinimumBiasHF1_OR_BptxAND',
    'L1_CastorHighJet_NotMinimumBiasHF2_AND_BptxAND',
    'L1_CastorHighJet_NotMinimumBiasHF2_OR_BptxAND',
    'L1_CastorHighJet_OR_MinimumBiasHF1_AND_BptxAND',
    'L1_CastorHighJet_OR_MinimumBiasHF2_AND_BptxAND',
    'L1_CastorMediumJet',
    'L1_CastorMediumJet_BptxAND',
    'L1_CastorMediumJet_MinimumBiasHF1_OR_BptxAND',
    'L1_CastorMediumJet_NotMinimumBiasHF2_AND_BptxAND',
    'L1_CastorMediumJet_NotMinimumBiasHF2_OR_BptxAND',
    'L1_CastorMediumJet_SingleEG5_MinimumBiasHF1_OR_BptxAND',
    'L1_CastorMediumJet_SingleMu0_MinimumBiasHF1_OR_BptxAND',
    'L1_CastorMuon',
    'L1_CastorMuon_BptxAND',
    'L1_Centrality_20_100_MinimumBiasHF1_AND_BptxAND',
    'L1_Centrality_30_100',
    'L1_Centrality_30_100_MinimumBiasHF1_AND_BptxAND',
    'L1_Centrality_50_100',
    'L1_Centrality_Saturation',
    'L1_DoubleEG10_BptxAND',
    'L1_DoubleEG2',
    'L1_DoubleEG2_BptxAND',
    'L1_DoubleEG2_NotMinimumBiasHF2_AND_BptxAND',
    'L1_DoubleEG2_NotMinimumBiasHF2_OR_BptxAND',
    'L1_DoubleEG5',
    'L1_DoubleEG5_BptxAND',
    'L1_DoubleEG5_NotMinimumBiasHF2_AND_BptxAND',
    'L1_DoubleEG5_NotMinimumBiasHF2_OR_BptxAND',
    'L1_DoubleEG8_BptxAND',
    'L1_DoubleJet16And12_MidEta2p7_BptxAND',
    'L1_DoubleJet16And12_MidEta2p7_Centrality_30_100_BptxAND',
    'L1_DoubleJet16And12_MidEta2p7_Centrality_50_100_BptxAND',
    'L1_DoubleJet16And8_MidEta2p7_BptxAND',
    'L1_DoubleJet16And8_MidEta2p7_Centrality_30_100_BptxAND',
    'L1_DoubleJet16And8_MidEta2p7_Centrality_50_100_BptxAND',
    'L1_DoubleJet20And12_MidEta2p7_BptxAND',
    'L1_DoubleJet20And12_MidEta2p7_Centrality_30_100_BptxAND',
    'L1_DoubleJet20And12_MidEta2p7_Centrality_50_100_BptxAND',
    'L1_DoubleJet20And8_MidEta2p7_BptxAND',
    'L1_DoubleJet20And8_MidEta2p7_Centrality_30_100_BptxAND',
    'L1_DoubleJet20And8_MidEta2p7_Centrality_50_100_BptxAND',
    'L1_DoubleJet28And16_MidEta2p7_BptxAND',
    'L1_DoubleJet28And16_MidEta2p7_Centrality_30_100_BptxAND',
    'L1_DoubleJet28And16_MidEta2p7_Centrality_50_100_BptxAND',
    'L1_DoubleMu0',
    'L1_DoubleMu0_BptxAND',
    'L1_DoubleMu0_Centrality_10_100_MinimumBiasHF1_AND_BptxAND',
    'L1_DoubleMu0_Centrality_30_100_MinimumBiasHF1_AND_BptxAND',
    'L1_DoubleMu0_Centrality_50_100_MinimumBiasHF1_AND_BptxAND',
    'L1_DoubleMu0_Mass_Min1',
    'L1_DoubleMu0_MinimumBiasHF1_AND_BptxAND',
    'L1_DoubleMu0_NotMinimumBiasHF2_AND_BptxAND',
    'L1_DoubleMu0_NotMinimumBiasHF2_OR_BptxAND',
    'L1_DoubleMu0_SQ',
    'L1_DoubleMu0_SQ_OS',
    'L1_DoubleMu10_BptxAND',
    'L1_DoubleMuOpen',
    'L1_DoubleMuOpen_BptxAND',
    'L1_DoubleMuOpen_Centrality_10_100_BptxAND',
    'L1_DoubleMuOpen_Centrality_30_100_BptxAND',
    'L1_DoubleMuOpen_Centrality_40_100_BptxAND',
    'L1_DoubleMuOpen_Centrality_50_100_BptxAND',
    'L1_DoubleMuOpen_MaxDr2p0_BptxAND',
    'L1_DoubleMuOpen_MaxDr2p0_OS_BptxAND',
    'L1_DoubleMuOpen_MaxDr3p5',
    'L1_DoubleMuOpen_MaxDr3p5_BptxAND',
    'L1_DoubleMuOpen_NotMinimumBiasHF2_AND_BptxAND',
    'L1_DoubleMuOpen_NotMinimumBiasHF2_OR_BptxAND',
    'L1_DoubleMuOpen_OS',
    'L1_DoubleMuOpen_OS_BptxAND',
    'L1_DoubleMuOpen_SS',
    'L1_DoubleMuOpen_SS_BptxAND',
    'L1_ETMHF100',
    'L1_ETMHF100_HTT60er',
    'L1_ETMHF120',
    'L1_ETMHF120_HTT60er',
    'L1_ETT10_ETTAsym50_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT10_ETTAsym55_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT10_ETTAsym60_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT10_ETTAsym65_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT10_ETTAsym70_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT10_ETTAsym80_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT1200',
    'L1_ETT1600',
    'L1_ETT2000',
    'L1_ETT35_NotETT80_BptxAND',
    'L1_ETT40_NotETT95_BptxAND',
    'L1_ETT45_NotETT110_BptxAND',
    'L1_ETT5',
    'L1_ETT50_ETTAsym40_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym40_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym50_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym50_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym55_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym60_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym60_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym65_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym70_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym70_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym80_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_ETTAsym80_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT50_NotETT120_BptxAND',
    'L1_ETT55_NotETT130_BptxAND',
    'L1_ETT5_BptxAND',
    'L1_ETT5_ETTAsym40_BptxAND',
    'L1_ETT5_ETTAsym40_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT5_ETTAsym40_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT5_ETTAsym50_BptxAND',
    'L1_ETT5_ETTAsym50_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT5_ETTAsym50_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT5_ETTAsym55_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT5_ETTAsym60_BptxAND',
    'L1_ETT5_ETTAsym60_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT5_ETTAsym60_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT5_ETTAsym65_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT5_ETTAsym70_BptxAND',
    'L1_ETT5_ETTAsym70_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT5_ETTAsym70_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT5_ETTAsym80_BptxAND',
    'L1_ETT5_ETTAsym80_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT5_ETTAsym80_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETT5_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT5_NotETT30_BptxAND',
    'L1_ETT5_NotMinimumBiasHF2_OR',
    'L1_ETT60_ETTAsym60_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT60_ETTAsym65_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT65_ETTAsym70_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT65_ETTAsym80_MinimumBiasHF2_OR_BptxAND',
    'L1_ETT8_ETTAsym50_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT8_ETTAsym55_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT8_ETTAsym60_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT8_ETTAsym65_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT8_ETTAsym70_MinimumBiasHF1_OR_BptxAND',
    'L1_ETT8_ETTAsym80_MinimumBiasHF1_OR_BptxAND',
    'L1_ETTAsym40',
    'L1_ETTAsym40_BptxAND',
    'L1_ETTAsym40_MinimumBiasHF1_OR_BptxAND',
    'L1_ETTAsym40_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETTAsym50',
    'L1_ETTAsym50_BptxAND',
    'L1_ETTAsym50_MinimumBiasHF1_OR_BptxAND',
    'L1_ETTAsym50_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETTAsym60',
    'L1_ETTAsym60_BptxAND',
    'L1_ETTAsym60_MinimumBiasHF1_OR_BptxAND',
    'L1_ETTAsym60_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETTAsym70',
    'L1_ETTAsym70_BptxAND',
    'L1_ETTAsym70_MinimumBiasHF1_OR_BptxAND',
    'L1_ETTAsym70_NotMinimumBiasHF2_OR_BptxAND',
    'L1_ETTAsym80',
    'L1_ETTAsym80_BptxAND',
    'L1_ETTAsym80_MinimumBiasHF1_OR_BptxAND',
    'L1_ETTAsym80_NotMinimumBiasHF2_OR_BptxAND',
    'L1_FirstBunchAfterTrain',
    'L1_FirstBunchBeforeTrain',
    'L1_FirstBunchInTrain',
    'L1_FirstCollisionInOrbit',
    'L1_FirstCollisionInOrbit_Centrality30_100_BptxAND',
    'L1_FirstCollisionInTrain',
    'L1_HCAL_LaserMon_Trig',
    'L1_HCAL_LaserMon_Veto',
    'L1_HTT120er',
    'L1_HTT200er',
    'L1_HTT280er',
    'L1_HTT360er',
    'L1_HTT450er',
    'L1_IsolatedBunch',
    'L1_LastBunchInTrain',
    'L1_LastCollisionInTrain',
    'L1_MinimumBiasHF0_AND_BptxAND',
    'L1_MinimumBiasHF0_OR_BptxAND',
    'L1_MinimumBiasHF1_AND',
    'L1_MinimumBiasHF1_AND_BptxAND',
    'L1_MinimumBiasHF1_AND_OR_ETT10_BptxAND',
    'L1_MinimumBiasHF1_OR',
    'L1_MinimumBiasHF1_OR_BptxAND',
    'L1_MinimumBiasHF1_XOR_BptxAND',
    'L1_MinimumBiasHF2_AND',
    'L1_MinimumBiasHF2_AND_BptxAND',
    'L1_MinimumBiasHF2_OR',
    'L1_MinimumBiasHF2_OR_BptxAND',
    'L1_NotBptxOR',
    'L1_NotETT100_MinimumBiasHF1_AND_BptxAND',
    'L1_NotETT110_MinimumBiasHF1_OR_BptxAND',
    'L1_NotETT110_MinimumBiasHF2_OR_BptxAND',
    'L1_NotETT150_MinimumBiasHF1_AND_BptxAND',
    'L1_NotETT150_MinimumBiasHF1_OR_BptxAND',
    'L1_NotETT150_MinimumBiasHF2_OR_BptxAND',
    'L1_NotETT200_MinimumBiasHF1_AND_BptxAND',
    'L1_NotETT20_MinimumBiasHF1_AND_BptxAND',
    'L1_NotETT20_MinimumBiasHF1_OR_BptxAND',
    'L1_NotETT20_MinimumBiasHF2_OR_BptxAND',
    'L1_NotETT80_MinimumBiasHF1_AND_BptxAND',
    'L1_NotETT80_MinimumBiasHF1_OR_BptxAND',
    'L1_NotETT80_MinimumBiasHF2_OR_BptxAND',
    'L1_NotETT95_MinimumBiasHF1_AND_BptxAND',
    'L1_NotETT95_MinimumBiasHF1_OR_BptxAND',
    'L1_NotETT95_MinimumBiasHF2_OR_BptxAND',
    'L1_NotMinimumBiasHF0_AND_BptxAND',
    'L1_NotMinimumBiasHF0_AND_BptxAND_TOTEM_1',
    'L1_NotMinimumBiasHF0_AND_BptxAND_TOTEM_2',
    'L1_NotMinimumBiasHF0_AND_BptxAND_TOTEM_4',
    'L1_NotMinimumBiasHF0_OR_BptxAND',
    'L1_NotMinimumBiasHF0_OR_BptxAND_TOTEM_1',
    'L1_NotMinimumBiasHF0_OR_BptxAND_TOTEM_2',
    'L1_NotMinimumBiasHF0_OR_BptxAND_TOTEM_4',
    'L1_NotMinimumBiasHF1_AND',
    'L1_NotMinimumBiasHF1_OR',
    'L1_NotMinimumBiasHF1_OR_BptxAND',
    'L1_NotMinimumBiasHF2_AND',
    'L1_NotMinimumBiasHF2_AND_BptxAND',
    'L1_NotMinimumBiasHF2_OR',
    'L1_NotMinimumBiasHF2_OR_BptxAND',
    'L1_SecondBunchInTrain',
    'L1_SecondLastBunchInTrain',
    'L1_SingleEG10er2p5',
    'L1_SingleEG12_BptxAND',
    'L1_SingleEG12_SingleJet28_MidEta2p7_BptxAND',
    'L1_SingleEG12_SingleJet28_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG12_SingleJet32_MidEta2p7_BptxAND',
    'L1_SingleEG12_SingleJet40_MidEta2p7_BptxAND',
    'L1_SingleEG12_SingleJet44_MidEta2p7_BptxAND',
    'L1_SingleEG12_SingleJet44_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG12_SingleJet56_MidEta2p7_BptxAND',
    'L1_SingleEG12_SingleJet56_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG12_SingleJet60_MidEta2p7_BptxAND',
    'L1_SingleEG12_SingleJet60_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG15_BptxAND',
    'L1_SingleEG15_Centrality_30_100_BptxAND',
    'L1_SingleEG15_Centrality_50_100_BptxAND',
    'L1_SingleEG15_SingleJet28_MidEta2p7_BptxAND',
    'L1_SingleEG15_SingleJet28_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG15_SingleJet44_MidEta2p7_BptxAND',
    'L1_SingleEG15_SingleJet44_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG15_SingleJet56_MidEta2p7_BptxAND',
    'L1_SingleEG15_SingleJet56_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG15_SingleJet60_MidEta2p7_BptxAND',
    'L1_SingleEG15_SingleJet60_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG15er2p5',
    'L1_SingleEG21_BptxAND',
    'L1_SingleEG21_Centrality_30_100_BptxAND',
    'L1_SingleEG21_Centrality_50_100_BptxAND',
    'L1_SingleEG26er2p5',
    'L1_SingleEG3',
    'L1_SingleEG30_BptxAND',
    'L1_SingleEG3_BptxAND',
    'L1_SingleEG3_Centrality_30_100_BptxAND',
    'L1_SingleEG3_Centrality_50_100_BptxAND',
    'L1_SingleEG3_NotMinimumBiasHF2_AND_BptxAND',
    'L1_SingleEG3_NotMinimumBiasHF2_OR_BptxAND',
    'L1_SingleEG5',
    'L1_SingleEG50',
    'L1_SingleEG5_BptxAND',
    'L1_SingleEG5_NotMinimumBiasHF2_AND_BptxAND',
    'L1_SingleEG5_NotMinimumBiasHF2_OR_BptxAND',
    'L1_SingleEG5_SingleJet28_MidEta2p7_BptxAND',
    'L1_SingleEG5_SingleJet32_MidEta2p7_BptxAND',
    'L1_SingleEG5_SingleJet40_MidEta2p7_BptxAND',
    'L1_SingleEG7_BptxAND',
    'L1_SingleEG7_Centrality_30_100_BptxAND',
    'L1_SingleEG7_Centrality_50_100_BptxAND',
    'L1_SingleEG7_SingleJet28_MidEta2p7_BptxAND',
    'L1_SingleEG7_SingleJet28_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG7_SingleJet32_MidEta2p7_BptxAND',
    'L1_SingleEG7_SingleJet40_MidEta2p7_BptxAND',
    'L1_SingleEG7_SingleJet44_MidEta2p7_BptxAND',
    'L1_SingleEG7_SingleJet44_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG7_SingleJet56_MidEta2p7_BptxAND',
    'L1_SingleEG7_SingleJet56_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG7_SingleJet60_MidEta2p7_BptxAND',
    'L1_SingleEG7_SingleJet60_MidEta2p7_MinDr0p4_BptxAND',
    'L1_SingleEG8er2p5',
    'L1_SingleIsoEG12_BptxAND',
    'L1_SingleIsoEG15_BptxAND',
    'L1_SingleIsoEG21_BptxAND',
    'L1_SingleIsoEG3_BptxAND',
    'L1_SingleIsoEG7_BptxAND',
    'L1_SingleJet120',
    'L1_SingleJet120_FWD3p0',
    'L1_SingleJet120er2p5',
    'L1_SingleJet16_BptxAND',
    'L1_SingleJet16_Centrality_30_100_BptxAND',
    'L1_SingleJet16_Centrality_50_100_BptxAND',
    'L1_SingleJet16_FWD_BptxAND',
    'L1_SingleJet16_FWD_Centrality_30_100_BptxAND',
    'L1_SingleJet16_FWD_Centrality_50_100_BptxAND',
    'L1_SingleJet180er2p5',
    'L1_SingleJet200',
    'L1_SingleJet24_BptxAND',
    'L1_SingleJet24_Centrality_30_100_BptxAND',
    'L1_SingleJet24_Centrality_50_100_BptxAND',
    'L1_SingleJet28_BptxAND',
    'L1_SingleJet28_Centrality_30_100_BptxAND',
    'L1_SingleJet28_Centrality_50_100_BptxAND',
    'L1_SingleJet28_FWD_BptxAND',
    'L1_SingleJet28_FWD_Centrality_30_100_BptxAND',
    'L1_SingleJet28_FWD_Centrality_50_100_BptxAND',
    'L1_SingleJet32_BptxAND',
    'L1_SingleJet32_Centrality_30_100_BptxAND',
    'L1_SingleJet32_Centrality_50_100_BptxAND',
    'L1_SingleJet35',
    'L1_SingleJet35_FWD3p0',
    'L1_SingleJet35er2p5',
    'L1_SingleJet36_BptxAND',
    'L1_SingleJet36_Centrality_30_100_BptxAND',
    'L1_SingleJet36_Centrality_50_100_BptxAND',
    'L1_SingleJet36_FWD_BptxAND',
    'L1_SingleJet36_FWD_Centrality_30_100_BptxAND',
    'L1_SingleJet36_FWD_Centrality_50_100_BptxAND',
    'L1_SingleJet40_BptxAND',
    'L1_SingleJet40_Centrality_30_100_BptxAND',
    'L1_SingleJet40_Centrality_50_100_BptxAND',
    'L1_SingleJet44_BptxAND',
    'L1_SingleJet44_Centrality_30_100_BptxAND',
    'L1_SingleJet44_Centrality_50_100_BptxAND',
    'L1_SingleJet44_FWD_BptxAND',
    'L1_SingleJet44_FWD_Centrality_30_100_BptxAND',
    'L1_SingleJet44_FWD_Centrality_50_100_BptxAND',
    'L1_SingleJet48_BptxAND',
    'L1_SingleJet48_Centrality_30_100_BptxAND',
    'L1_SingleJet48_Centrality_50_100_BptxAND',
    'L1_SingleJet56_BptxAND',
    'L1_SingleJet56_Centrality_30_100_BptxAND',
    'L1_SingleJet56_Centrality_50_100_BptxAND',
    'L1_SingleJet56_FWD_BptxAND',
    'L1_SingleJet56_FWD_Centrality_30_100_BptxAND',
    'L1_SingleJet56_FWD_Centrality_50_100_BptxAND',
    'L1_SingleJet60',
    'L1_SingleJet60_BptxAND',
    'L1_SingleJet60_Centrality_30_100_BptxAND',
    'L1_SingleJet60_Centrality_50_100_BptxAND',
    'L1_SingleJet60_FWD3p0',
    'L1_SingleJet60er2p5',
    'L1_SingleJet64_BptxAND',
    'L1_SingleJet64_Centrality_30_100_BptxAND',
    'L1_SingleJet64_Centrality_50_100_BptxAND',
    'L1_SingleJet64_FWD_BptxAND',
    'L1_SingleJet64_FWD_Centrality_30_100_BptxAND',
    'L1_SingleJet64_FWD_Centrality_50_100_BptxAND',
    'L1_SingleJet72_BptxAND',
    'L1_SingleJet8',
    'L1_SingleJet80_BptxAND',
    'L1_SingleJet8_BptxAND',
    'L1_SingleJet8_Centrality_30_100_BptxAND',
    'L1_SingleJet8_Centrality_50_100_BptxAND',
    'L1_SingleJet8_FWD_BptxAND',
    'L1_SingleJet8_FWD_Centrality_30_100_BptxAND',
    'L1_SingleJet8_FWD_Centrality_50_100_BptxAND',
    'L1_SingleJet90',
    'L1_SingleJet90_FWD3p0',
    'L1_SingleJet90er2p5',
    'L1_SingleMu0',
    'L1_SingleMu0_BMTF',
    'L1_SingleMu0_BptxAND',
    'L1_SingleMu0_DQ',
    'L1_SingleMu0_EMTF',
    'L1_SingleMu0_NotMinimumBiasHF2_AND_BptxAND',
    'L1_SingleMu0_NotMinimumBiasHF2_OR_BptxAND',
    'L1_SingleMu0_OMTF',
    'L1_SingleMu12',
    'L1_SingleMu12_BptxAND',
    'L1_SingleMu12_DQ_BMTF',
    'L1_SingleMu12_DQ_EMTF',
    'L1_SingleMu12_DQ_OMTF',
    'L1_SingleMu12_MinimumBiasHF1_AND_BptxAND',
    'L1_SingleMu12_SingleEG7_BptxAND',
    'L1_SingleMu16',
    'L1_SingleMu16_BptxAND',
    'L1_SingleMu16_MinimumBiasHF1_AND_BptxAND',
    'L1_SingleMu22',
    'L1_SingleMu22_BMTF',
    'L1_SingleMu22_EMTF',
    'L1_SingleMu22_OMTF',
    'L1_SingleMu3',
    'L1_SingleMu3Open_BptxAND',
    'L1_SingleMu3_BptxAND',
    'L1_SingleMu3_Centrality_70_100_BptxAND',
    'L1_SingleMu3_Centrality_80_100_BptxAND',
    'L1_SingleMu3_MinimumBiasHF1_AND_BptxAND',
    'L1_SingleMu3_NotMinimumBiasHF2_OR_BptxAND',
    'L1_SingleMu3_SingleEG12_BptxAND',
    'L1_SingleMu3_SingleEG15_BptxAND',
    'L1_SingleMu3_SingleEG20_BptxAND',
    'L1_SingleMu3_SingleEG30_BptxAND',
    'L1_SingleMu3_SingleJet28_MidEta2p7_BptxAND',
    'L1_SingleMu3_SingleJet32_MidEta2p7_BptxAND',
    'L1_SingleMu3_SingleJet40_MidEta2p7_BptxAND',
    'L1_SingleMu5',
    'L1_SingleMu5_BptxAND',
    'L1_SingleMu5_MinimumBiasHF1_AND_BptxAND',
    'L1_SingleMu5_SingleEG10_BptxAND',
    'L1_SingleMu5_SingleEG12_BptxAND',
    'L1_SingleMu5_SingleEG15_BptxAND',
    'L1_SingleMu5_SingleEG20_BptxAND',
    'L1_SingleMu7',
    'L1_SingleMu7_BptxAND',
    'L1_SingleMu7_DQ',
    'L1_SingleMu7_MinimumBiasHF1_AND_BptxAND',
    'L1_SingleMu7_SingleEG10_BptxAND',
    'L1_SingleMu7_SingleEG12_BptxAND',
    'L1_SingleMu7_SingleEG15_BptxAND',
    'L1_SingleMu7_SingleEG7_BptxAND',
    'L1_SingleMuCosmics',
    'L1_SingleMuCosmics_BMTF',
    'L1_SingleMuCosmics_EMTF',
    'L1_SingleMuCosmics_OMTF',
    'L1_SingleMuOpen',
    'L1_SingleMuOpen_BptxAND',
    'L1_SingleMuOpen_Centrality_70_100_BptxAND',
    'L1_SingleMuOpen_Centrality_70_100_MinimumBiasHF1_AND_BptxAND',
    'L1_SingleMuOpen_Centrality_80_100_BptxAND',
    'L1_SingleMuOpen_Centrality_80_100_MinimumBiasHF1_AND_BptxAND',
    'L1_SingleMuOpen_NotMinimumBiasHF2_AND_BptxAND',
    'L1_SingleMuOpen_NotMinimumBiasHF2_OR_BptxAND',
    'L1_SingleMuOpen_SingleEG15_BptxAND',
    'L1_SingleMuOpen_SingleJet28_MidEta2p7_BptxAND',
    'L1_SingleMuOpen_SingleJet44_MidEta2p7_BptxAND',
    'L1_SingleMuOpen_SingleJet56_MidEta2p7_BptxAND',
    'L1_SingleMuOpen_SingleJet64_MidEta2p7_BptxAND',
    'L1_TOTEM_1',
    'L1_TOTEM_2',
    'L1_TOTEM_3',
    'L1_TOTEM_4',
    'L1_UnpairedBunchBptxMinus',
    'L1_UnpairedBunchBptxPlus',
    'L1_ZDCM',
    'L1_ZDCM_BptxAND',
    'L1_ZDCM_ZDCP_BptxAND',
    'L1_ZDCP',
    'L1_ZDCP_BptxAND',
    'L1_ZDC_AND_OR_MinimumBiasHF1_AND_BptxAND',
    'L1_ZDC_AND_OR_MinimumBiasHF1_OR_BptxAND',
    'L1_ZDC_AND_OR_MinimumBiasHF2_AND_BptxAND',
    'L1_ZDC_AND_OR_MinimumBiasHF2_OR_BptxAND',
    'L1_ZDC_OR_OR_MinimumBiasHF1_OR_BptxAND',
    'L1_ZDC_OR_OR_MinimumBiasHF2_OR_BptxAND',
    'L1_ZeroBias',
    'L1_ZeroBias_copy',
    ])
