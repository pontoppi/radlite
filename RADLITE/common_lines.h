#ifdef INCLUDE_LINES
      common/lines/Aud,Aud_const,Bud,Bdu,Jbarline,linefreq,dek,gratio
      doubleprecision Aud(SZ_NLINES),Aud_const(SZ_NLINES)
      doubleprecision Bud(SZ_NLINES),Bdu(SZ_NLINES)
      doubleprecision Jbarline(SZ_NLINES)
      doubleprecision linefreq(SZ_NLINES),dek(SZ_NLINES)
      doubleprecision gratio(SZ_NLINES)
      common/ilines/lev_up,lev_down,nlines,nlin_orig
      integer lev_up(SZ_NLINES),lev_down(SZ_NLINES)
      integer nlines,nlin_orig
c
      common/colltrans/Cud,Cdu,Kud,TempTrans,gratio_ct,dek_ct
      doubleprecision TempTrans(SZ_NCTTEMP)
      doubleprecision Cud(SZ_NCOLTRANS)
      doubleprecision Cdu(SZ_NCOLTRANS)
      doubleprecision Kud(SZ_NCOLTRANS,SZ_NCTTEMP)
      doubleprecision gratio_ct(SZ_NCOLTRANS)
      doubleprecision dek_ct(SZ_NCOLTRANS)
      common/icolltrans/ncoltrans,ncoltr_orig,ncttemp,
     %                  lev_ct_up,lev_ct_down
      integer ncoltrans,ncoltr_orig,ncttemp
      integer lev_ct_up(SZ_NCOLTRANS)
      integer lev_ct_down(SZ_NCOLTRANS)
c
      common/levels/npopul,gdeg,enerlevel,enerlevcm
      doubleprecision npopul(SZ_NLEVELS),gdeg(SZ_NLEVELS)
      doubleprecision enerlevel(SZ_NLEVELS),enerlevcm(SZ_NLEVELS)
      common/ilevels/nlevels,nlev_orig,ifop
      integer nlevels,nlev_orig,ifop
c
      common/linedum/Jbarplanck
      doubleprecision Jbarplanck(SZ_NLINES)
c
      common/lineprofiles/line_phi,line_nu0,line_wgt,line_rangewidth,
     %                    line_fix_linewidth,line_dnu,line_nu,line_phi2,
     %                    aksmin,aksmax
      doubleprecision line_phi(SZ_LINEPROFILE,SZ_NLINES)
      doubleprecision line_nu0(SZ_NLINES)
      doubleprecision line_dnu(SZ_LINEPROFILE,SZ_NLINES)
      doubleprecision line_nu(SZ_LINEPROFILE,SZ_NLINES)
      doubleprecision line_wgt(SZ_LINEPROFILE,SZ_NLINES)
      doubleprecision line_rangewidth(SZ_NLINES)
      doubleprecision line_fix_linewidth(SZ_NLINES)
      doubleprecision line_phi2(SZ_LINEPROFILE,SZ_NLINES)
      doubleprecision aksmin,aksmax
      common/ilineprofiles/line_symm,line_nfreq,line_profiletype
      integer line_symm(SZ_NLINES),line_nfreq(SZ_NLINES)
      integer line_profiletype(SZ_NLINES)
c
      common/linklinetrans/line_trifreq,line_iline,line_ilifreq
      integer line_trifreq(SZ_LINEPROFILE,SZ_NLINES)
      integer line_iline(FRSIZE_FREQ),line_ilifreq(FRSIZE_FREQ)
c
      common/linefield/line_jbar_field,line_level_popul,
     %     Aud_field,Bud_field,Bdu_field,Cud_field,Cdu_field,
     %     line_approx_lambdabar,line_source_func,
     %     line_jbar_field_prev
      doubleprecision Aud_field(1:SZ_NLINES,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision Bud_field(1:SZ_NLINES,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision Bdu_field(1:SZ_NLINES,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision Cud_field(1:SZ_NCOLTRANS,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision Cdu_field(1:SZ_NCOLTRANS,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision line_jbar_field(1:SZ_NLINES,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision line_jbar_field_prev(1:SZ_NLINES,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision line_level_popul(1:SZ_NLEVELS,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision line_approx_lambdabar(1:SZ_NLINES,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
      doubleprecision line_source_func(1:SZ_NLINES,
     %                 0:FRSIZE_Y_SMALL,0:FRSIZE_X)
c
      common/linewarn/line_warn_lineout,line_warn_cr_extrapol,
     %                line_warn_negpop,line_warn_badsampling,
     %                line_warn_preptrans_done,line_warn_velostep,
     %                line_warn_opacity_neg,line_warn_neg_opacity,
     %                line_pop_lastsave,line_warn_dust_done,
     %                line_jbar_lastsave,line_warn_umass
      integer line_warn_lineout,line_warn_cr_extrapol
      integer line_warn_negpop,line_warn_badsampling
      integer line_warn_preptrans_done,line_warn_velostep
      integer line_warn_opacity_neg,line_warn_neg_opacity
      integer line_pop_lastsave,line_warn_dust_done
      integer line_jbar_lastsave,line_warn_umass
c
      common/lineaprlambda/lvalamb_om_nu
      doubleprecision lvalamb_om_nu(SZ_LINEPROFILE_SMALL,
     %      1:FRSIZE_PHI_SMALL,-FRSIZE_MU_HALF:FRSIZE_MU_HALF)
c
c     If LINE_VELOCITIES is specified, the line transfer can be done using
c     all the sophistication possible: systematic velocities, local
c     abundances, local line width etc. If this is not set, then the
c     line transfer is done in a primitive way, with only a global line
c     profile, global abundance, global line width etc.
c
#ifdef LINE_VELOCITIES
c
c     If LINE_LOCAL_WIDTH is set, then the line profile width can be
c     different at different locations. This function a(R,Theta) should
c     then be specified by the setup routine. If LINE_LOCAL_WIDTH is
c     not set, then the line width is taken from the line.inp file.
c
#ifdef LINE_LOCAL_WIDTH
c
c     If LINE_PROFILE_ARRAY is set, then the local line profiles are stored
c     in advance in a big array. This has the advantage of speed, but the
c     disadvantage of large memory requirements.
c
      common/localwidths/locprof_linewidth
      doubleprecision locprof_linewidth(1:SZ_NLINES,
     %                     1:FRSIZE_Y_SMALL,1:FRSIZE_X)
#ifdef LINE_PROFILE_ARRAY
      common/localprofiles/locprof_phi,locprof_phi2
      doubleprecision locprof_phi(SZ_LINEPROFILE,SZ_NLINES,
     %                     1:FRSIZE_Y_SMALL,1:FRSIZE_X)
      doubleprecision locprof_phi2(SZ_LINEPROFILE,SZ_NLINES,
     %                     1:FRSIZE_Y_SMALL,1:FRSIZE_X)
#endif
#endif
c
c     If LINE_LOCAL_ABUNDANCE is set, then the abundance is assumed
c     to be specified in the array locabun_abund_mol(), which should
c     be filled by the setup routine. If not set, the abundance
c     is assumed to be global, and specified in the line.inp file.
c
#ifdef LINE_LOCAL_ABUNDANCE
      common/localabundance/locabun_abund_mol,locabun_abund_coll
      doubleprecision locabun_abund_mol(1:FRSIZE_Y_SMALL,
     %                                  1:FRSIZE_X)
      doubleprecision locabun_abund_coll(1:FRSIZE_Y_SMALL,
     %                                  1:FRSIZE_X)
      common/ilocalabundance/locabun_done
      integer locabun_done
#endif
#endif
c
c     It is advised always to define LINE_NG_ACCEL, if line transfer
c     is enabled. This enables the crucial Ng acceleration!
c
c     Some stuff relevant to dust in line transfer
c
      common/linedust/line_dust_src,line_dust_alp
#ifndef LINE_FREQDEP_DUST
      doubleprecision line_dust_src(1,SZ_NLINES,
     %                     1:FRSIZE_Y_SMALL,1:FRSIZE_X)
      doubleprecision line_dust_alp(1,SZ_NLINES,
     %                     1:FRSIZE_Y_SMALL,1:FRSIZE_X)
#else /* LINE_FREQDEP_DUST */
      doubleprecision line_dust_src(SZ_LINEPROFILE,SZ_NLINES,
     %                     1:FRSIZE_Y_SMALL,1:FRSIZE_X)
      doubleprecision line_dust_alp(SZ_LINEPROFILE,SZ_NLINES,
     %                     1:FRSIZE_Y_SMALL,1:FRSIZE_X)
#endif /* LINE_FREQDEP_DUST */
c
c     NEW 27.03.06: for combined dust and lines
c     
      common/ilinedust/iline_freq_nr,maserflag
      integer iline_freq_nr
      integer maserflag
c
      common/linepump/line_starmi,line_extmi
      doubleprecision line_starmi(SZ_NLINES),line_extmi(SZ_NLINES)
c
#endif  /* IFDEF INCLUDE_LINES */
c
      common/molecules/umass_av,umass_molec,umass_colpartner,
     %                          abund_molec,abund_colpartner
      doubleprecision umass_av,umass_molec,umass_colpartner
      doubleprecision abund_molec,abund_colpartner
      common/molecname/molec_name,molec_filename
      character*200 molec_name,molec_filename
      common/imolecname/molec_name_len
      integer molec_name_len
c
