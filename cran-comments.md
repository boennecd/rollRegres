## Test environments
* Platform: x86_64-w64-mingw32/x64 (64-bit)
  Running under: Windows >= 8 x64 (build 9200)
  R version 3.5.0
* Ubuntu 14.04.5 LTS (on travis-ci with codename: trusty)
  R version 3.5.0
* Local Ubuntu 18.04 (64-bit) with R devel and with clang 6.0.0 with ASAN and 
  UBSAN checks
* The following rhub platforms:
  debian-gcc-devel
  debian-gcc-patched
  debian-gcc-release
  fedora-clang-devel
  fedora-gcc-devel
  solaris-x86-patched
  linux-x86_64-rocker-gcc-san
* win-builder (devel and release)

## R CMD check results
No errors and only the note about `New submission`

## Resubmission
This is a resubmission. In this version I have:

* Changed the description. I have added a reference as Swetlana Herbrandt
  requested. 
* Fix these two bugs: 
  https://github.com/boennecd/rollRegres/issues/1 and 
  https://github.com/boennecd/rollRegres/issues/2 
