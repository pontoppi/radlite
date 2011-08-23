c     --------------------------------------------------------------
c                THE GLOBAL COMMONS FOR THE DUST TEMPERATURES
c                              AND OPACITIES
c     --------------------------------------------------------------
#ifdef INCLUDE_DUST
#ifndef INCLUDE_DUST_ARRAYS
#define INCLUDE_DUST_ARRAYS
#endif
#endif
c
#ifdef INCLUDE_DUST_ARRAYS
      common/idustopac/dust_done_read,dust_opacity_tempdep
      integer dust_done_read,dust_opacity_tempdep
      common/dustopac/dust_kappawgt_abs,dust_kappawgt_scat,
     %                dust_freq_wgt,dust_tmin,dust_tmax,
     %                dust_dtmin,dust_dtmax,dust_temprange_high,
     %                dust_temprange_low,dust_mgrain
      doubleprecision dust_kappawgt_abs(1:FRSIZE_FREQ,1:DUST_TRANGE_MAX,1:DUST_SIZE_MAX,1:DUST_SPECIES_MAX)
      doubleprecision dust_kappawgt_scat(1:FRSIZE_FREQ,1:DUST_TRANGE_MAX,
     %                    1:DUST_SIZE_MAX,1:DUST_SPECIES_MAX)
      doubleprecision dust_freq_wgt(1:FRSIZE_FREQ)
      doubleprecision dust_temprange_low(1:DUST_TRANGE_MAX,
     %                    1:DUST_SPECIES_MAX)
      doubleprecision dust_temprange_high(1:DUST_TRANGE_MAX,
     %                    1:DUST_SPECIES_MAX)
      doubleprecision dust_tmin(1:DUST_SPECIES_MAX)
      doubleprecision dust_tmax(1:DUST_SPECIES_MAX)
      doubleprecision dust_dtmin(1:DUST_SPECIES_MAX)
      doubleprecision dust_dtmax(1:DUST_SPECIES_MAX)
      doubleprecision dust_mgrain(1:DUST_SPECIES_MAX)
c
#     ifdef DUST_RADICAL
      common/dustdenstemp/dust_temp,dust_rho,ngbk_dust_temp
      doubleprecision dust_rho(1:DUST_SPECIES_MAX,
     %                    0:FRSIZE_Y,0:FRSIZE_X) 
      doubleprecision dust_temp(1:DUST_SIZE_MAX,1:DUST_SPECIES_MAX,
     %                    0:FRSIZE_Y,0:FRSIZE_X) 
#     ifndef NO_DUST_NG
      doubleprecision ngbk_dust_temp(1:DUST_SIZE_MAX,1:DUST_SPECIES_MAX,
     %                    0:FRSIZE_Y,0:FRSIZE_X,4)
#     else
      doubleprecision ngbk_dust_temp
#     endif
#     endif
c
      common/idust/dust_nr_species,dust_nr_size,
     %             dust_frwgt_read,dust_warn_few_freqs,
     %             dust_nr_temp,dust_warn_zero_temp,
     %             dust_quantum
      integer dust_nr_species,dust_warn_zero_temp
      integer dust_frwgt_read,dust_warn_few_freqs
      integer dust_nr_size(1:DUST_SPECIES_MAX)
      integer dust_nr_temp(1:DUST_SPECIES_MAX)
      integer dust_quantum(1:DUST_SPECIES_MAX)
c
#ifdef DUST_TRIDAG_ALI
      common/dusttriali/dust_diaglamb_a,dust_diaglamb_a
      doubleprecision dust_diaglamb_a(1:FRSIZE_Y_SMALL,1:FRISZE_X)
      doubleprecision dust_diaglamb_s(1:FRSIZE_Y_SMALL,1:FRISZE_X)
c
      common/dustmeankap/jmeankap_a,jmeankap_s,pmeankap_a,pmeankap_s,
     %                   imeankap_a,imeankap_s,imeankap_denom
      doubleprecision jmeankap_a(1:FRSIZE_Y_SMALL,1:FRISZE_X)
      doubleprecision jmeankap_s(1:FRSIZE_Y_SMALL,1:FRISZE_X)
      doubleprecision pmeankap_a(1:FRSIZE_Y_SMALL,1:FRISZE_X)
      doubleprecision pmeankap_s(1:FRSIZE_Y_SMALL,1:FRISZE_X)
      doubleprecision imeankap_a(1:FRSIZE_PHI_SMALL,
     %     -FRSIZE_MU_HALF:FRSIZE_MU_HALF,1:FRSIZE_Y_SMALL,1:FRISZE_X)
      doubleprecision imeankap_s(1:FRSIZE_PHI_SMALL,
     %     -FRSIZE_MU_HALF:FRSIZE_MU_HALF,1:FRSIZE_Y_SMALL,1:FRISZE_X)
      doubleprecision imeankap_denom(1:FRSIZE_PHI_SMALL,
     %     -FRSIZE_MU_HALF:FRSIZE_MU_HALF,1:FRSIZE_Y_SMALL,1:FRISZE_X)
#endif
c
#ifdef DUST_ALI
      common/dustali/dust_tempold,dust_meanintold
      doubleprecision dust_tempold(1:FRSIZE_Y_SMALL,1:FRSIZE_X)
      doubleprecision dust_meanintold(1:FRSIZE_FREQ,
     %                   1:FRSIZE_Y_SMALL,1:FRSIZE_X)
#endif
c
#endif
c
      common/idustsetup/dust_setup_nrspecies,dust_setup_nrsizes
      integer dust_setup_nrspecies
      integer dust_setup_nrsizes(1:DUST_SPECIES_MAX)
c
#ifdef DUST_PRIVATE_EXTERNALS
      common/freqgrid/freq_nu
      doubleprecision freq_nu(FRSIZE_FREQ)
      common/ifreqgrid/freq_nr
      integer freq_nr
#endif
c
