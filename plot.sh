rivet-mkhtml \
  -m 'njl' \
  --single \
  --errs \
  fused.yoda:'Title=fused (detector uncertainty)':'BandComponentEnv=btagrate,ctagrate,jes':'Variations=none' \
  fused1.yoda:'Title=fused1' \
  -o test.plots
  # 'fused1.yoda:Title=fused ($\mu_R$ up):DefaultWeight=MUR2_MUF1_PDF261000:Variations=none' \
  # 'fused2.yoda:Title=fused ($\mu_R$ down):DefaultWeight=MUR0.5_MUF1_PDF261000:Variations=none' \
  # 'direct.yoda:Title=direct:Variations=none' \
  # 'fragmentation.yoda:Title=fragmentation:Variations=none' \
  # 'fused.yoda:Title=fused:Variations=btagrate,ctagrate,jes,MUR2_MUF2_PDF261000,MUR0.5_MUF0.5_PDF261000' \
  # 'direct.yoda:Title=direct:Variations=btagrate,ctagrate,jes,MUR2_MUF2_PDF261000,MUR0.5_MUF0.5_PDF261000' \
  # 'fragmentation.yoda:Title=fragmentation:Variations=btagrate,ctagrate,jes,MUR2_MUF2_PDF261000,MUR0.5_MUF0.5_PDF261000'
