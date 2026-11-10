# Chain A ... HF-RPA & HF-RRPA sources in dependency order (dispatcher last) ...
SRC_RPA = src/rpa/nj_hf_rpa.f90 \
          src/rpa/nj_hf_rpa_export.f90 \
          src/rpa/nj_hf_rpa_external.f90 \
          src/rpa/nj_hf_rpa_driver.f90 \
          src/rpa/nj_hf_rrpa.f90 \
          src/rpa/nj_hf_rrpa_export.f90 \
          src/rpa/nj_hf_rrpa_driver.f90 \
          src/rpa/nj_dispatch_rpa.f90
