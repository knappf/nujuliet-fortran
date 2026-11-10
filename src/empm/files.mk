# Chain C ... EMPM (restricted-space F-matrix, TDA, 2- & 3-phonon EMPM) sources in dependency order (dispatcher last) ...
# Uses chain A (nj_hf_rpa, nj_hf_rpa_export) ... build it with the stub dispatcher when chain A is stubbed out ...
SRC_EMPM = src/empm/nj_empm_space.f90 \
           src/empm/nj_empm_fmat.f90 \
           src/empm/nj_empm_tda.f90 \
           src/empm/nj_empm_export.f90 \
           src/empm/nj_empm_phonons.f90 \
           src/empm/nj_empm_dens.f90 \
           src/empm/nj_empm_eqm2.f90 \
           src/empm/nj_empm_full.f90 \
           src/empm/nj_empm_eqm3.f90 \
           src/empm/nj_empm_cm.f90 \
           src/empm/nj_empm_driver.f90 \
           src/empm/nj_dispatch_empm.f90
