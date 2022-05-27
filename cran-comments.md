## Test environments

*local: Windows 11 
*r_hub:
*check_mac_release()


## R CMD check

1. check()
── R CMD check results ─────────────────── ASV 1.0 ────
Duration: 5m 4.9s

❯ checking top-level files ... NOTE
  Non-standard files/directories found at top level:
    'Readme.md' 'cran-comments.md'

0 errors ✔ | 0 warnings ✔ | 1 note ✖

2. check_rhub()
── ASV 1.0: NOTE

  Build ID:   ASV_1.0.tar.gz-e88acfc2cea64757be95f27fdd0c116f
  Platform:   Windows Server 2022, R-devel, 64 bit
  Submitted:  11m 42.2s ago
  Build time: 11m 40.2s

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Yaushiro Omori <omori.yasuhiro@gmail.com>'
  
  New submission

❯ checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
               user system elapsed
  ReportMCMC  81.39   0.20   81.74
  sv_mcmc     78.55   0.05   78.61
  asv_apf     53.99   1.89   42.14
  asv_pf      27.73   0.69   24.94
  sv_pf       27.54   0.67   25.12
  asv_mcmc   180.20   0.19  192.03
  sv_apf      27.16   0.89   25.00

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

── ASV 1.0: IN-PROGRESS

  Build ID:   ASV_1.0.tar.gz-cbed15d6fa164a3da3ae6feee15e6dd9
  Platform:   Ubuntu Linux 20.04.1 LTS, R-release, GCC
  Submitted:  11m 42.3s ago


── ASV 1.0: IN-PROGRESS

  Build ID:   ASV_1.0.tar.gz-fb6000821d4445a8885484997c0a9bea
  Platform:   Fedora Linux, R-devel, clang, gfortran
  Submitted:  11m 42.3s ago


── ASV 1.0: OK

  Build ID:   ASV_1.0.tar.gz-f4eaffc9daf14b0a94757509acc20826
  Platform:   Debian Linux, R-devel, GCC ASAN/UBSAN
  Submitted:  11m 42.3s ago
  Build time: 11m 5s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

3.check_mac_release()
✔ [2022-05-27] Check <https://mac.R-project.org/macbuilder/results/1653612385-812893a32072715b/> for the results in 5-10 mins (~09:56).


## revdepcheck results

There are currently no downstream dependencies for this package

